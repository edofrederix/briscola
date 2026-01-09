#include "FadlunDirichletImmersedBoundaryCondition.H"
#include "immersedBoundary.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

template<class Type, class MeshType>
FadlunDirichletImmersedBoundaryCondition<Type,MeshType>::
FadlunDirichletImmersedBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const immersedBoundary<MeshType>& ib
)
:
    immersedBoundaryCondition<Type,MeshType>(field, ib, &ib.wallAdjMask()),
    boundaryValues_(this->read("value"))
{
    // Check shape overlap

    if (this->ib_.shapeOverlap())
        this->shapeOverlapWarning();

    // Check for closely packed shapes

    const meshField<label,MeshType>& mask = this->forcingMask();

    forAllCells(mask,l,d,i,j,k)
    if (mask(l,d,i,j,k))
    {
        const labelVector ijk(i,j,k);

        for (int f = 0; f < 3; f++)
        {
            const labelVector fo = units[f];

            if
            (
                this->ib_.mask()[l][d](ijk + fo)
             && this->ib_.mask()[l][d](ijk - fo)
            )
            {
                WarningInFunction
                    << "Wall adjacent cell " << ijk
                    << " at level " << l << ", direction " << d
                    << " has an immersed boundary on both sides."
                    << endl;

                break;
            }
        }
    }
}

// Destructor

template<class Type, class MeshType>
FadlunDirichletImmersedBoundaryCondition<Type,MeshType>::
~FadlunDirichletImmersedBoundaryCondition()
{}

template<class Type, class MeshType>
void FadlunDirichletImmersedBoundaryCondition<Type,MeshType>::evaluate
(
    const label l,
    const label d
)
{
    const scalar omega = this->omega_;

    meshDirection<Type,MeshType>& x = this->field_[l][d];

    const meshDirection<label,MeshType>& mask = this->forcingMask()[l][d];
    const meshDirection<faceScalar,MeshType>& y = this->ib_.wallDistAdj()[l][d];

    forAllCells(x,i,j,k)
    {
        const labelVector ijk(i,j,k);

        if (mask(ijk))
        {
            scalar ximax = 0;

            // Loop over face number directions
            for (int f = 0; f < 6; f++)
            {
                const labelVector fo = faceOffsets[f];

                if (y(ijk)[f] > ximax)
                {
                    ximax = y(ijk)[f];

                    const scalar xic = 1.0 - ximax;
                    const scalar xinb = 1.0 + xic;
                    const scalar w = xic/xinb;

                    Type value =
                        boundaryValues_[d]/xinb + w*x(ijk - fo);

                    x(ijk) = (1.0 - omega)*x(ijk) + omega*value;
                }
            }
        }
    }
}

}

}

}
