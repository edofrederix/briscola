#include "immersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Return colocated mass source term to be added to the Poisson equation
// given a staggered velocity field. The source term corrects for unphysical
// fluxes across the immersed boundary that can occur in methods that use
// velocity forcing.

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

    forAll(fvMsh.IB<staggered>(), ib)
    {
        forAllCells(ghostMask,l,d,i,j,k)
        {
            if (fvMsh.IB<staggered>()[ib].ghostMask()(l,d,i,j,k))
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

    forAllCells(source[0][0],i,j,k)
    {
        source(0,0,i,j,k) -=
            ghostMask(0,0,i,j,k) * field(0,0,i,j,k) * fa(0,0,i,j,k).left();

        source(0,0,i,j,k) +=
            ghostMask(0,0,i+1,j,k) * field(0,0,i+1,j,k) * fa(0,0,i,j,k).right();

        source(0,0,i,j,k) -=
            ghostMask(0,1,i,j,k) * field(0,1,i,j,k) * fa(0,0,i,j,k).bottom();

        source(0,0,i,j,k) +=
            ghostMask(0,1,i,j+1,k) * field(0,1,i,j+1,k) * fa(0,0,i,j,k).top();

        source(0,0,i,j,k) -=
            ghostMask(0,2,i,j,k) * field(0,2,i,j,k) * fa(0,0,i,j,k).aft();

        source(0,0,i,j,k) +=
            ghostMask(0,2,i,j,k+1) * field(0,2,i,j,k+1) * fa(0,0,i,j,k).fore();

        source(0,0,i,j,k) /= cv(0,0,i,j,k);
    }

    return tSource;
}

}

}

}