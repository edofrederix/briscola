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

template<class ViscosityModel>
void twoPhaseMultiVof<ViscosityModel>::initAlphas()
{
    alphas_.clear();

    label index = 0;

    if
    (
        this->fvMsh().db().template foundObject<colocatedScalarField>
        (
            "alpha.0"
        )
    )
    {
        while
        (
            this->fvMsh().db().template foundObject<colocatedScalarField>
            (
                "alpha."
                +Foam::name(index)
            )
        )
        {
            colocatedScalarField& alpha = this->fvMsh().db().template
                lookupObjectRef<colocatedScalarField>
                (
                    "alpha."
                    +Foam::name(index)
                );

            // Create tmp copy of alpha field
            tmp<colocatedScalarField> tAlpha(new colocatedScalarField(alpha));

            // Release and delete old alpha field
            alpha.regIOobject::release();
            delete &alpha;

            // Create new vofField with tmp alpha field
            alphas_.append
            (
                new vofField
                (
                    "alpha."+Foam::name(index),
                    tAlpha(),
                    this->dict(),
                    false
                )
            );

            index++;
        }
    }
    else
    {
        alphas_.append
        (
            new vofField
            (
                "alpha.0",
                this->alpha_,
                this->dict()
            )
        );
    }

}

template<class ViscosityModel>
void twoPhaseMultiVof<ViscosityModel>::addAlphaField()
{
    alphas_.append
    (
        new vofField
        (
            "alpha."+Foam::name(alphas_.size()),
            this->alpha_,
            this->dict(),
            true
        )
    );

    #ifdef FULLDEBUG

    Info << "Alpha field added" << endl;

    #endif
}

template<class ViscosityModel>
void twoPhaseMultiVof<ViscosityModel>::computeGlobaltagAlgorithm()
{
    tags_ = Zero;

    forAllCells(tags_, i,j,k)
    {
        forAll(alphas_, a)
        {
            if (alphas_[a].tag()(i,j,k))
            {
                tags_(i,j,k) = alphas_[a].tag()(i,j,k);

                break;
            }
        }
    }
}

template<class ViscosityModel>
void twoPhaseMultiVof<ViscosityModel>::computeGlobalN()
{
    N_ = Zero;
    scalar removedVol = Zero;

    forAll(alphas_, a)
    {
        N_ += alphas_[a].tag().n();
        removedVol += alphas_[a].tag().removedVol();
    }

    Pstream::gather(removedVol, sumOp<scalar>());

    Info << "Total number of particles: " << N_
         << " in " << alphas_.size() << " alpha fields. "
         << "Total volume of small particles removed: "
         << removedVol << endl;
}


template<class ViscosityModel>
void twoPhaseMultiVof<ViscosityModel>::setPhi()
{
    phi_.clear();

    forAll(alphas_, a)
    {
        for (int i = 0; i < alphas_[a].tag().n(); i++)
        {
            phi_.append(a);
        }
    }
}


template<class ViscosityModel>
void twoPhaseMultiVof<ViscosityModel>::overwriteTagsWithIDs()
{
    int Np = 0;

    forAll(alphas_, a)
    {
        forAllCells(alphas_[a].tag(),i,j,k)
        {
            if (alphas_[a].tag()(i,j,k))
            {
                alphas_[a].tag()(i,j,k) += Np;
            }
        }

        Np += alphas_[a].tag().n();
    }
}

template<class ViscosityModel>
void twoPhaseMultiVof<ViscosityModel>::setConnectivityMatrix()
{
    connectivityMatrix_.setSize(N_);

    for (int i = 0; i < connectivityMatrix_.m(); i++)
        for (int j = 0; j < connectivityMatrix_.n(); j++)
            connectivityMatrix_(i,j) = false;

    forAllCells(tags_, i,j,k)
    {
        labelList tagNums(0);

        for (int di = -1; di <= 1; di++)
        {
            for (int dj = -1; dj <= 1; dj++)
            {
                for (int dk = -1; dk <= 1; dk++)
                {
                    label tagNum = tags_(i + di, j + dj, k + dk);

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

                    connectivityMatrix_(first-1, second-1) = true;
                    connectivityMatrix_(second-1, first-1) = true;

                }
            }
        }
    }

    for (int i = 0; i < connectivityMatrix_.m(); i++)
        for (int j = 0; j < connectivityMatrix_.n(); j++)
            reduce(connectivityMatrix_(i,j), orOp<bool>());
}

template<class ViscosityModel>
void twoPhaseMultiVof<ViscosityModel>::moveFields()
{
    for (int i = 0; i < connectivityMatrix_.m()-1; i++)
    {
        for (int j = i+1; j < connectivityMatrix_.n(); j++)
        {
            if (connectivityMatrix_(i,j) && phi_[i] == phi_[j])
            {
                // Create ineligible interface array. This array identifies
                // which interfaces are ineligible for hosting the particle.

                List<bool> ineligibleAlphas(alphas_.size(), false);

                // The current interface number is ineligible, by definition

                ineligibleAlphas[phi_[i]] = true;

                // Fields that have particles that are connected to the current
                // one are ineligible too

                for (int m = 0; m < N_; m++)
                {
                    label phi = phi_[m];

                    if (phi != phi_[i] && connectivityMatrix_(i,m))
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

                Info << "Moving particle " << i+1 << " from field "
                     << phi_[i] << " to field " << targetAlpha << endl;

                #endif

                // Move the particle from the source interface to the target interface

                forAllCells(alphas_[phi_[i]],x,y,z)
                {
                    if (alphas_[phi_[i]].tag()(x,y,z) == i+1)
                    {
                        alphas_[targetAlpha].tag()(x,y,z) = i+1;
                        alphas_[phi_[i]].tag()(x,y,z) = 0;

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

                for (int m = 0; m < connectivityMatrix_.m(); m++)
                {
                    connectivityMatrix_[i][m] = false;
                    connectivityMatrix_[m][i] = false;
                }
            }
        }
    }
}

template<class ViscosityModel>
twoPhaseMultiVof<ViscosityModel>::twoPhaseMultiVof
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    ViscosityModel(fvMsh, dict),
    alphas_(),
    tags_
    (
        "tags",
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
    N_(-1),
    phi_(),
    connectivityMatrix_()
{
    tags_ = Zero;
}

template<class ViscosityModel>
twoPhaseMultiVof<ViscosityModel>::twoPhaseMultiVof
(
    const twoPhaseMultiVof& tpm
)
:
    ViscosityModel(tpm),
    alphas_(tpm.alphas_),
    tags_(tpm.tags_),
    colors_(tpm.colors_),
    N_(tpm.N_),
    phi_(tpm.phi_),
    connectivityMatrix_(tpm.connectivityMatrix_)
{
    forAll(alphas_, a)
    {
        alphas_[a].normal().correct();
        alphas_[a].surfaceTension().correct();
    }

    ViscosityModel::correct();
}

template<class ViscosityModel>
tmp<colocatedScalarFaceField>
twoPhaseMultiVof<ViscosityModel>::flux()
{
    tmp<colocatedScalarFaceField> tSurfaceTensionFlux
    (
        new colocatedScalarFaceField("surfaceTensionFlux", this->fvMsh_)
    );

    colocatedScalarFaceField& flux = tSurfaceTensionFlux.ref();

    flux = Zero;

    forAll(alphas_, a)
    {
        flux += static_cast<colocatedScalarFaceField&>
            (
                alphas_[a].surfaceTension()
            );
    }

    return tSurfaceTensionFlux;
}

template<class ViscosityModel>
void twoPhaseMultiVof<ViscosityModel>::correct()
{
    if (alphas_.empty())
    {
        initAlphas();
    }

    // Tag alpha fields
    forAll(alphas_, a)
    {
        alphas_[a].tag().correct();
    }

    // Coalescence suppression
    computeGlobalN();
    setPhi();
    overwriteTagsWithIDs();
    computeGlobaltagAlgorithm();
    setConnectivityMatrix();
    moveFields();

    // Set colors field
    colors_ = tags_;

    forAllCells(colors_,i,j,k)
    {
        if (colors_(i,j,k))
        {
            colors_(i,j,k) = phi_[colors_(i,j,k)-1] + 1;
        }
    }

    // Correct multiple marker fields
    forAll(alphas_, a)
    {
        alphas_[a].vf().solve(this->coloFaceFlux()());
    }

    // Compute summed alpha field

    this->alpha_ = Zero;

    forAll(alphas_, a)
    {
        this->alpha_ += static_cast<colocatedScalarField&>(alphas_[a]);
    }

    // Set bounds for summed alpha field
    this->alpha_.correctAlpha();

    // Correct the mesh-specific face volume fraction
    ViscosityModel::correctFaceAlpha();

    // Correct the base model
    ViscosityModel::correct();

    // Correct surface tensions
    forAll(alphas_, a)
    {
        alphas_[a].surfaceTension().correct();
    }
}

}

}

}
