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

tmp<colocatedScalarField> ibmCorr
(
    tmp<colocatedScalarField> tColoDivU,
    const staggeredScalarField& field
)
{
    const fvMesh& fvMsh = field.fvMsh();

    Switch massSource = false;

    forAll(field.immersedBoundaryConditions(), i)
    {
        const word IBname = fvMsh.ibs<staggered>()[i].name();

        if
        (
            field.subDict("boundaryConditions").subDict(IBname)
                .lookupOrDefault<Switch>("massSource", false)
        )
        {
            massSource = true;
        }
    }

    colocatedScalarField& coloDivU = tColoDivU.ref();

    if (massSource)
    {
        meshField<label,staggered> ghostMask
        (
            "massSourceGhostMask",
            fvMsh
        );

        ghostMask = Zero;

        forAll(fvMsh.ibs<staggered>(), ib)
        {
            forAllCells(ghostMask,l,d,i,j,k)
            {
                if (fvMsh.ibs<staggered>()[ib].ghostMask()(l,d,i,j,k))
                {
                    ghostMask(l,d,i,j,k) = 1;
                }
            }
        }

        const colocatedFaceScalarField& fa =
            fvMsh.metrics<colocated>().faceAreas();

        const colocatedScalarField& cv =
            fvMsh.metrics<colocated>().cellVolumes();

        forAllCells(coloDivU, i, j, k)
        {
            const labelVector ijk(i,j,k);

            for (int d = 0; d < 3; d++)
            {
                const labelVector nei(upperNei(ijk,d));

                coloDivU(ijk) +=
                    (
                        ghostMask(ijk)*field(ijk)*fa(ijk)[d*2]
                      - ghostMask(nei)*field(nei)*fa(ijk)[d*2+1]
                    ) / cv(ijk);
            }
        }
    }

    return tColoDivU;
}

}

}

}