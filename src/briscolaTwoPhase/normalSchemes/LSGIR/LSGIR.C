#include "LSGIR.H"
#include "addToRunTimeSelectionTable.H"
#include "vof.H"
#include "exSchemes.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(LSGIR, 0);
addToRunTimeSelectionTable(normalScheme, LSGIR, dictionary);

LSGIR::LSGIR
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const colocatedScalarField& alpha
)
:
    normalScheme(fvMsh, dict, alpha),
    beta_(dict.lookupOrDefault<scalar>("beta", 1.5))
{}

LSGIR::LSGIR(const LSGIR& s)
:
    normalScheme(s),
    beta_(s.beta_)
{}

LSGIR::~LSGIR()
{}

void LSGIR::correct()
{
    colocatedVectorField& n = *this;

    // Reconstruct the interface normal using the version of LSGIR presented in
    // Lopez & Hernández (2022)

    const meshField<vector,colocated>& cc =
        fvMsh_.template metrics<colocated>().cellCenters();

    const faceLabel I(n.I());

    forAllCells(n, i, j, k)
    {
        const labelVector ijk(i,j,k);

        if
        (
            (alpha_(ijk) > vof::threshold)
         && (alpha_(ijk) < 1 - vof::threshold)
        )
        {
            // Loop over 3^3-1 neighboring cells and compute A and b matrices

            vectorList A(26, Zero);
            scalarList b(26, Zero);

            label c = 0;
            labelVector o;
            for (o.x() = -1; o.x() < 2; o.x()++)
            for (o.y() = -1; o.y() < 2; o.y()++)
            for (o.z() = -1; o.z() < 2; o.z()++)
            if (o != zeroXYZ)
            {
                // Check if we are in a neighboring cell that is internal or is
                // a parallel/periodic ghost cell

                const bool internal =
                    cmptMin(ijk + o - I.lower()) > -1
                 && cmptMax(ijk + o - I.upper()) <  0;

                if
                (
                    internal
                 || fvMsh_.msh().boundaries().parallelMask()(o + unitXYZ)
                )
                {
                    const scalar weight =
                        1.0/Foam::pow(mag(cc(ijk) - cc(ijk+o)), beta_);

                    A[c] = weight*(cc(ijk+o) - cc(ijk));
                    b[c] = weight*(alpha_(ijk+o) - alpha_(ijk));

                    c++;
                }
            }

            A.resize(c);
            b.resize(c);

            // Compute the matrix product A^TA and matrix-vector product A^Tb

            tensor ATA = Zero;
            vector ATb = Zero;

            forAll(A, l)
                for (label ii = 0; ii < 3; ii++)
                    for (label jj = 0; jj < 3; jj++)
                        ATA(ii,jj) += A[l][ii]*A[l][jj];

            forAll(b, l)
                for (label ii = 0; ii < 3; ii++)
                    ATb[ii] += A[l][ii]*b[l];

            // Solve stabilized linear system

            for (label ii = 0; ii < 3; ii++)
                if (trimPrecision(ATA(ii,ii)) == 0)
                    ATA(ii,ii) = SMALL;

            n(ijk) = ATA.inv() & ATb;

            const scalar magn = Foam::mag(n(ijk));

            if (magn > 0.0)
                n(ijk) /= magn;
        }
        else
        {
            n(ijk) = Zero;
        }
    }

    n.correctBoundaryConditions();
}

}

}

}
