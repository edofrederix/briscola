#include "immersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Return colocated mass source term to be added to the Poisson equation given a
// staggered velocity field. The source term corrects for non-physical fluxes
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
        const word ibName = fvMsh.immersedBoundaries<staggered>()[i].name();

        if
        (
            field.subDict("boundaryConditions").subDict(ibName)
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

        forAll(fvMsh.immersedBoundaries<staggered>(), ib)
        {
            const immersedBoundary<staggered>& b =
                fvMsh.immersedBoundaries<staggered>()[ib];

            forAllCells(ghostMask,l,d,i,j,k)
            {
                if (b.ghostMask()(l,d,i,j,k))
                {
                    ghostMask(l,d,i,j,k) = 1;
                }
            }
        }

        const colocatedScalarFaceField& fa =
            fvMsh.metrics<colocated>().faceAreas();

        const colocatedScalarField& icv =
            fvMsh.metrics<colocated>().inverseCellVolumes();

        coloDivU = Zero;

        forAllFaces(fa, fd, i, j, k)
        {
            const labelVector ijk(i,j,k);
            const labelVector nei(lowerNeighbor(i,j,k,fd));

            const scalar value =
                ghostMask(fd,ijk)*field(fd,ijk)*fa[fd](ijk);

            coloDivU(ijk) += value;
            coloDivU(nei) -= value;
        }

        coloDivU *= icv;
    }

    return tColoDivU;
}

}

}

}