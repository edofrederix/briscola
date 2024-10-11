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
    const meshField<Type,MeshType>& mshField,
    const immersedBoundary<MeshType>& ib
)
:
    immersedBoundaryCondition<Type,MeshType>
    (
        mshField,
        ib,
        ib.wallAdjMask()
    ),
    boundaryValues_(this->dict().lookup("values"))
{
    // Check shape overlap
    if (this->IB_.shapeOverlap())
    {
        WarningInFunction
            << "Overlapping shapes identified."
            << " This may cause issues with Fadlun IBM." << endl;
    }

    // Check for closely packed shapes
    forAllCells(this->IB_.wallAdjMask(),l,d,i,j,k)
    {
        const labelVector ijk(i,j,k);

        if (this->IB_.wallAdjMask()(l,d,i,j,k))
        {
            for (int dir = 0; dir < 3; dir++)
            {
                const labelVector fo = faceOffsets[2*dir];

                if
                (
                       this->IB_.mask()[l][d](ijk+fo)
                    && this->IB_.mask()[l][d](ijk-fo)
                )
                {
                    WarningInFunction
                        << "Wall adjacent cell "
                        << vector(i,j,k) << " at (l,d) = "<< l << ", "
                        << d << " has immersed boundary on both sides."
                        << " This may cause issues with Fadlun IBM."
                        << endl;

                    break;
                }
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
void FadlunDirichletImmersedBoundaryCondition<Type,MeshType>::
correctJacobiPoints
(
    meshLevel<Type,MeshType>& x
) const
{
    scalar omega = this->omega_;

    label l = x.levelNum();

    if (l == 0)
    {
        forAllCells(x,d,i,j,k)
        {
            if (this->forcingPoints_(l,d,i,j,k))
            {
                scalar ximax = 0;

                // Loop over face number directions
                for (int dir = 0; dir < 6; dir++)
                {
                    if (this->IB_.wallDistAdj()(l,d,i,j,k)[dir] > ximax)
                    {
                        ximax = this->IB_.wallDistAdj()(l,d,i,j,k)[dir];
                        const scalar xic = 1.0
                            - this->IB_.wallDistAdj()(l,d,i,j,k)[dir];
                        const scalar xinb = 1.0 + xic;
                        const scalar w = xic/xinb;

                        labelVector ijk(i,j,k);

                        Type forcingValue = boundaryValues_[d]/xinb
                            + w*x[d](ijk-faceOffsets[dir]);

                        x(d,i,j,k) = (1.0 - omega) * x(d,i,j,k)
                            + omega * forcingValue;
                    }
                }
            }
        }
    }
}

}

}

}
