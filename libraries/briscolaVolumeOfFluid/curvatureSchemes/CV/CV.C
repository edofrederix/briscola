#include "CV.H"
#include "addToRunTimeSelectionTable.H"
#include "vof.H"
#include "exSchemes.H"
#include "rectilinearMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(CV, 0);
addToRunTimeSelectionTable(curvatureScheme, CV, dictionary);

CV::CV(const vof& vf, const dictionary& dict)
:
    curvatureScheme(vf, dict)
{}

CV::CV(const CV& s)
:
    curvatureScheme(s)
{}

CV::~CV()
{}

tmp<colocatedScalarField> CV::operator()(const colocatedVectorField& n)
{
    tmp<colocatedScalarField> tconvolutedAlpha
    (
        new colocatedScalarField
        (
            "alpha",
            vf_.fvMsh()
        )
    );

    colocatedScalarField& convolutedAlpha = tconvolutedAlpha.ref();

    const colocatedScalarField& alpha = vf_.alpha();

    const meshField<vector,colocated>& centers =
        vf_.fvMsh().template metrics<colocated>().cellCenters();

    scalar h = Foam::mag(centers(1,1,1) - centers(0,0,0));
    reduce(h, maxOp<scalar>());

    double kernel[3][3] = {
        {1.0/16.0, 2.0/16.0, 1.0/16.0},
        {2.0/16.0, 4.0/16.0, 2.0/16.0},
        {1.0/16.0, 2.0/16.0, 1.0/16.0}
    };

    forAllCells(alpha, i, j, k)
    {
        convolutedAlpha(i,j,k) = 0;

        for (int aux1 = -1; aux1 <= 1; aux1++)
        {
            for (int aux2 = -1; aux2 <= 1; aux2++)
            {
                    convolutedAlpha(i,j,k) +=
                        kernel[aux1 + 1][aux2 + 1]
                        *alpha(i+aux1,j+aux2,k);
            }
        }
    }

    tmp<colocatedVectorField> tn
    (
        new colocatedVectorField(ex::grad(tconvolutedAlpha()))
    );

    colocatedVectorField& normal = tn.ref();

    forAllCells(alpha, i, j, k)
    {
        scalar S = Foam::mag(normal(i,j,k));

        if (S > 1e-8)
            normal(i,j,k) /= S;
        else
            normal(i,j,k) = Zero;
    }

    const colocatedFaceScalarField nFlux
    (
        ex::faceFlux(normal)
    );

    tmp<colocatedScalarField> tk2
    (
        new colocatedScalarField(ex::div(nFlux))
    );

    colocatedScalarField& kc2 = tk2.ref();

    tmp<colocatedScalarField> tk
    (
        new colocatedScalarField
        (
            "curvature",
            vf_.fvMsh()
        )
    );

    colocatedScalarField& kc = tk.ref();

    forAllCells(alpha, i, j, k)
    {
        if
        (
            (alpha(i,j,k) > vof::threshold)
            && (alpha(i,j,k) < 1 - vof::threshold)
        )
        {
            int count = 0;
            kc(i,j,k) = 0;

            for (int aux1 = -1; aux1 <= 1; aux1++)
            {
                for (int aux2 = -1; aux2 <= 1; aux2++)
                {
                    if
                    (
                        (alpha(i+aux1,j+aux2,k) > vof::threshold)
                        && (alpha(i+aux1,j+aux2,k) < 1 - vof::threshold)
                    )
                    {
                        count++;
                        kc(i,j,k) += kc2(i+aux1,j+aux2,k);
                    }
                }
            }

            kc(i,j,k) /= double(count);
        }
        else
            kc(i,j,k) = 0;
    }

    return tk2;
}

}

}

}
