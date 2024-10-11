#include "immersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Return colocated mass source term to be added to the Poisson equation given a
// staggered velocity field. The source term corrects for unphysical fluxes
// across the immersed boundary that can occur in methods that use velocity
// forcing.

tmp<colocatedScalarField> IBMMassSource
(
    const staggeredScalarField& field
)
{
    const fvMesh& fvMsh = field.fvMsh();

    meshField<label,staggered> ghostMask
    (
        "massSourceGhostMask",
        fvMsh
    );

    ghostMask = Zero;

    forAll(fvMsh.IBs<staggered>(), ib)
    {
        forAllCells(ghostMask,l,d,i,j,k)
        {
            if (fvMsh.IBs<staggered>()[ib].ghostMask()(l,d,i,j,k))
            {
                ghostMask(l,d,i,j,k) = 1;
            }
        }
    }

    tmp<colocatedScalarField> tSource
    (
        new colocatedScalarField
        (
            "IBMSource",
            field.fvMsh()
        )
    );

    colocatedScalarField& source = tSource.ref();

    source = Zero;

    const colocatedFaceScalarField& fa =
        fvMsh.metrics<colocated>().faceAreas();

    const colocatedScalarField& cv =
        fvMsh.metrics<colocated>().cellVolumes();

    forAllCells(source, i, j, k)
    {
        const labelVector ijk(i,j,k);

        for (int d = 0; d < 3; d++)
        {
            const labelVector nei(upperNei(ijk,d));

            source(ijk) -=
                ghostMask(ijk)*field(ijk)*fa(ijk)[d*2]
              - ghostMask(nei)*field(nei)*fa(ijk)[d*2+1];
        }

        source(ijk) /= cv(ijk);
    }

    return tSource;
}

}

}

}