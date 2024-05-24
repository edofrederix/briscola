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
FadlunDirichletImmersedBoundaryCondition<Type,MeshType>
::FadlunDirichletImmersedBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const immersedBoundary<MeshType>& ib
)
:
    immersedBoundaryCondition<Type,MeshType>(mshField,ib),
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
FadlunDirichletImmersedBoundaryCondition<Type,MeshType>
::~FadlunDirichletImmersedBoundaryCondition()
{}

template<class Type, class MeshType>
void FadlunDirichletImmersedBoundaryCondition<Type,MeshType>
::correctLinearSystem
(
    linearSystem<stencil,Type,MeshType>& ls
)
{
    forAllCells(ls.b(),l,d,i,j,k)
    {
        if (this->IB_.wallAdjMask()(l,d,i,j,k))
        {
            // Set sources to 0 in IB
            ls.b()(l,d,i,j,k) = Zero;
        }
    }

    forAllCells(ls.A(),l,d,i,j,k)
    {
        if (this->IB_.wallAdjMask()(l,d,i,j,k))
        {
            scalar ximax = 0;
            // Loop over face number directions
            for (int dir = 0; dir < 6; dir++)
            {
                const label oppositeDir =
                    faceNumber(-faceOffsets[dir]);

                if (this->IB_.wallDistAdj()(l,d,i,j,k)[dir] > ximax)
                {
                    ximax = this->IB_.wallDistAdj()(l,d,i,j,k)[dir];
                    const scalar xic = 1.0
                        - this->IB_.wallDistAdj()(l,d,i,j,k)[dir];
                    const scalar xinb = 1.0 + xic;
                    const scalar w = xic/xinb;

                    ls.A()(l,d,i,j,k) = Zero;
                    ls.A()(l,d,i,j,k).center() = 1.0;
                    ls.A()(l,d,i,j,k)[oppositeDir+1] = -w;

                    if (l==0)
                    {
                        ls.b()(l,d,i,j,k) = boundaryValues_[d]/xinb;
                    }
                }
            }
        }
    }
}

}

}

}
