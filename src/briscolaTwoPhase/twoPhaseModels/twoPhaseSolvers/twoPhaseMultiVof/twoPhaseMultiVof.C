#include "twoPhaseMultiVof.H"

#include "addToRunTimeSelectionTable.H"
#include "exSchemes.H"
#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
void twoPhaseMultiVof<BaseModel>::initAlphas()
{
    alphas_.clear();

    alphas_.append
    (
        new colocatedScalarField
        (
            "alpha0",
            this->alpha_,
            true,
            true
        )
    );

    alphas_[0].setRestrictionScheme("volumeWeighted");

    normalSchemes_.clear();

    normalSchemes_.append
    (
        normalScheme::New
        (
            this->fvMsh(),
            this->dict().subDict("normalScheme"),
            alphas_[0]
        )
    );

    surfaceTensionSchemes_.clear();

    surfaceTensionSchemes_.append
    (
        surfaceTensionScheme::New
        (
            this->fvMsh(),
            this->dict().subDict("surfaceTensionScheme"),
            normalSchemes_[0],
            alphas_[0]
        )
    );

    vofs_.clear();

    vofs_.append
    (
        vof::New
        (
            this->fvMsh(),
            this->dict().subDict("vof"),
            normalSchemes_[0],
            alphas_[0]
        )
    );

    CCL_.clear();

    CCL_.append
    (
        CCL::New
        (
            this->fvMsh(),
            this->dict().found("CCL")
            ? this->dict().subDict("CCL")
            : dictionary::null,
            alphas_[0]
        )
    );
}

template<class BaseModel>
void twoPhaseMultiVof<BaseModel>::addAlphaField()
{
    alphas_.append
    (
        new colocatedScalarField
        (
            "alpha"+Foam::name(alphas_.size()),
            this->alpha_,
            true,
            true
        )
    );

    alphas_.last() = Zero;

    alphas_.last().setRestrictionScheme("volumeWeighted");

    normalSchemes_.append
    (
        normalScheme::New
        (
            this->fvMsh(),
            this->dict().subDict("normalScheme"),
            alphas_.last()
        )
    );

    surfaceTensionSchemes_.append
    (
        surfaceTensionScheme::New
        (
            this->fvMsh(),
            this->dict().subDict("surfaceTensionScheme"),
            normalSchemes_.last(),
            alphas_.last()
        )
    );

    vofs_.append
    (
        vof::New
        (
            this->fvMsh(),
            this->dict().subDict("vof"),
            normalSchemes_.last(),
            alphas_.last()
        )
    );

    CCL_.append
    (
        CCL::New
        (
            this->fvMsh(),
            this->dict().found("CCL")
            ? this->dict().subDict("CCL")
            : dictionary::null,
            alphas_.last()
        )
    );

    #ifdef FULLDEBUG

    Info << "Alpha field added" << endl;

    #endif
}

template<class BaseModel>
void twoPhaseMultiVof<BaseModel>::computeGlobalCCL()
{
    globalCCL_ = Zero;

    forAllCells(globalCCL_, i,j,k)
    {
        forAll(CCL_,c)
        {
            if (CCL_[c](i,j,k))
            {
                globalCCL_(i,j,k) = CCL_[c](i,j,k);

                break;
            }
        }
    }
}

template<class BaseModel>
void twoPhaseMultiVof<BaseModel>::computeGlobalN()
{
    globalN_ = Zero;

    forAll(CCL_, c)
    {
        globalN_ += CCL_[c].n();
    }

    Info << "Total number of connected components: " << globalN_
         << " in: " << CCL_.size() << " alpha fields" << endl;
}


template<class BaseModel>
void twoPhaseMultiVof<BaseModel>::setPhi()
{
    phi_.clear();

    forAll(alphas_, a)
    {
        for (int i = 0; i < CCL_[a].n(); i++)
        {
            phi_.append(a);
        }
    }
}


template<class BaseModel>
void twoPhaseMultiVof<BaseModel>::overwriteTagsWithIDs()
{
    int Np = 0;

    forAll(CCL_,c)
    {
        forAllCells(CCL_[c],i,j,k)
        {
            if (CCL_[c](i,j,k))
            {
                CCL_[c](i,j,k) += Np;
            }
        }

        Np += CCL_[c].n();
    }
}

template<class BaseModel>
void twoPhaseMultiVof<BaseModel>::setConnectivityMatrix()
{
    connectivityMatrix_.setSize(globalN_);

    forAll(connectivityMatrix_, i)
    {
        connectivityMatrix_[i].setSize(globalN_);

        forAll(connectivityMatrix_[i], j)
        {
            connectivityMatrix_[i][j] = false;
        }
    }

    forAllCells(globalCCL_, i,j,k)
    {
        labelList tagNums(0);

        for (int di = -1; di <= 1; di++)
        {
            for (int dj = -1; dj <= 1; dj++)
            {
                for (int dk = -1; dk <= 1; dk++)
                {
                    label tagNum = globalCCL_(i + di, j + dj, k + dk);

                    if (tagNum)
                    {
                        tagNums.append(tagNum);
                    }
                }
            }
        }

        forAll(tagNums, m)
        {
            forAll(tagNums, n)
            {
                if (tagNums[m] != tagNums[n])
                {
                    label first = tagNums[m];
                    label second = tagNums[n];

                    connectivityMatrix_[first-1][second-1] = true;
                    connectivityMatrix_[second-1][first-1] = true;

                }
            }
        }
    }

    forAll(connectivityMatrix_, i)
    {
        forAll(connectivityMatrix_[i], j)
        {
            reduce(connectivityMatrix_[i][j], orOp<bool>());
        }
    }
}

template<class BaseModel>
void twoPhaseMultiVof<BaseModel>::moveFields()
{
    for (int i = 0; i < connectivityMatrix_.size()-1; i++)
    {
        for (int j = i+1; j < connectivityMatrix_[i].size(); j++)
        {
            if (connectivityMatrix_[i][j] && phi_[i] == phi_[j])
            {
                // Create ineligible interface array. This array identifies
                // which interfaces are ineligible for hosting the particle.

                List<bool> ineligibleAlphas(alphas_.size(), false);

                // The current interface number is ineligible, by definition

                ineligibleAlphas[phi_[i]] = true;

                // Fields that have particles that are connected to the current
                // one are ineligible too

                for (int m = 0; m < globalN_; m++)
                {
                    label phi = phi_[m];

                    if (phi != phi_[i] && connectivityMatrix_[i][m])
                    {
                        ineligibleAlphas[phi] = true;
                    }
                }

                // Communicate the ineligible array by a collective 'or' operation

                forAll(ineligibleAlphas, a)
                {
                    reduce(ineligibleAlphas[a], orOp<bool>());
                }

                // Find eligible interface

                label targetAlpha = -1;

                forAll(ineligibleAlphas, a)
                {
                    if (!ineligibleAlphas[a])
                    {
                        targetAlpha = a;
                    }
                }

                // Create and append new interface if there is no eligible one.
                // The target interface becomes this newly created one.

                if (targetAlpha == -1)
                {
                    addAlphaField();

                    targetAlpha = alphas_.size() - 1;
                }

                #ifdef FULLDEBUG

                Info << "Moving bubble " << i+1 << " from field "
                     << phi_[i] << " to field " << targetAlpha << endl;

                #endif

                // Move the particle from the source interface to the target interface

                forAllCells(alphas_[phi_[i]],x,y,z)
                {
                    if (CCL_[phi_[i]](x,y,z) == i+1)
                    {
                        CCL_[targetAlpha](x,y,z) = i+1;
                        CCL_[phi_[i]](x,y,z) = 0;

                        alphas_[targetAlpha](x,y,z) = Foam::min
                            (
                                alphas_[targetAlpha](x,y,z)
                              + alphas_[phi_[i]](x,y,z),
                                1.0
                            );

                        alphas_[phi_[i]](x,y,z) = 0.0;
                    }
                }

                // Correct boundary conditions of both modified alpha fields

                alphas_[phi_[i]].correctBoundaryConditions();
                alphas_[targetAlpha].correctBoundaryConditions();

                // Update interface number
                phi_[i] = targetAlpha;

                // Update connectivity matrix. By definition, the particle is
                // moved in such a way that it is no longer connected to any
                // other particle, so the i's row and column in the connectivity
                // matrix should be false.

                forAll(connectivityMatrix_, m)
                {
                    connectivityMatrix_[i][m] = false;
                    connectivityMatrix_[m][i] = false;
                }
            }
        }
    }
}

template<class BaseModel>
twoPhaseMultiVof<BaseModel>::twoPhaseMultiVof
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    BaseModel(fvMsh, dict),
    alphas_(),
    normalSchemes_(),
    surfaceTensionSchemes_(),
    vofs_(),
    CCL_(),
    globalCCL_
    (
        "m",
        fvMsh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true,
        false
    ),
    colors_
    (
        "colors",
        fvMsh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true,
        false
    ),
    globalN_(-1),
    phi_(),
    connectivityMatrix_()
{
    globalCCL_ = Zero;
}

template<class BaseModel>
twoPhaseMultiVof<BaseModel>::twoPhaseMultiVof
(
    const twoPhaseMultiVof& tpm
)
:
    BaseModel(tpm),
    alphas_(tpm.alphas_),
    normalSchemes_(tpm.normalSchemes_),
    surfaceTensionSchemes_(tpm.surfaceTensionSchemes_),
    vofs_(tpm.vofs_),
    CCL_(tpm.CCL_),
    globalCCL_(tpm.globalCCL_),
    colors_(tpm.colors_),
    globalN_(tpm.globalN_),
    phi_(tpm.phi_),
    connectivityMatrix_(tpm.connectivityMatrix_)
{
    forAll(normalSchemes_,i)
    {
        normalSchemes_[i].correct();
    }

    forAll(surfaceTensionSchemes_,i)
    {
        surfaceTensionSchemes_[i].correct();
    }

    BaseModel::correctMixture();
}

template<class BaseModel>
twoPhaseMultiVof<BaseModel>::~twoPhaseMultiVof()
{}

template<class BaseModel>
tmp<colocatedFaceScalarField>
twoPhaseMultiVof<BaseModel>::surfaceTensionFlux()
{
    tmp<colocatedFaceScalarField> tSurfaceTensionFlux
    (
        new colocatedFaceScalarField("surfaceTensionFlux", this->fvMsh_)
    );

    colocatedFaceScalarField& flux = tSurfaceTensionFlux.ref();

    #ifdef NO_BLOCK_ZERO_INIT
    flux = Zero;
    #endif

    forAll(surfaceTensionSchemes_, i)
    {
        flux += static_cast<colocatedFaceScalarField&>
            (
                surfaceTensionSchemes_[i]
            );
    }

    return tSurfaceTensionFlux;
}

template<class BaseModel>
void twoPhaseMultiVof<BaseModel>::correct()
{
    if (alphas_.empty())
    {
        initAlphas();
    }

    // Tag alpha fields
    forAll(CCL_,i)
    {
        CCL_[i].correct();
    }

    // Coalescence suppression
    computeGlobalN();
    setPhi();
    overwriteTagsWithIDs();
    computeGlobalCCL();
    setConnectivityMatrix();
    moveFields();

    // Set colors field
    colors_ = globalCCL_;

    forAllCells(colors_,i,j,k)
    {
        if (colors_(i,j,k))
        {
            colors_(i,j,k) = phi_[colors_(i,j,k)-1] + 1;
        }
    }

    // Correct multiple marker fields
    forAll(vofs_,i)
    {
        vofs_[i].solve(this->coloFaceFlux()());
    }

    // Compute summed alpha field

    this->alpha_ = Zero;

    forAll(alphas_,i)
    {
        this->alpha_ += alphas_[i];
    }

    // Set bounds for summed alpha field

    this->alpha_.correctBoundaryConditions();

    scalarBlock& alpha = this->alpha_.B();

    forAllBlockLinear(alpha, i)
        alpha(i) =
            Foam::min
            (
                Foam::max
                (
                    round(0.5*alpha(i)/vof::threshold)*2.0*vof::threshold,
                    0.0
                ),
                1.0
            );

    // Correct surface tensions
    forAll(surfaceTensionSchemes_,i)
    {
        surfaceTensionSchemes_[i].correct();
    }

    // Correct interface normals
    forAll(normalSchemes_,i)
    {
        normalSchemes_[i].correct();
    }

    BaseModel::correctMixture();
}

}

}

}
