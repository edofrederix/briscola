#include "Fadlun.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

template<class Type, class MeshType>
Fadlun<Type,MeshType>::Fadlun
(
    dictionary& dict,
    const fvMesh& fvMsh
)
:
    immersedBoundaryMethod<Type,MeshType>(dict,fvMsh,false),
    wallAdjMask_
    (
        "wallAdjMask",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        true
    ),
    wallDist_
    (
        "wallDist",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        true
    )
{
    // Cell centers
    const meshField<vector,MeshType>& CC =
        fvMsh.metrics<MeshType>().cellCenters();

    // Set IB mask fields
    forAllCells(wallAdjMask_,l,d,i,j,k)
    {
        const labelVector ijk(i,j,k);

        wallAdjMask_(l,d,i,j,k) = 0;

        if (!this->isInside(CC(l,d,i,j,k)))
        {
            const vector c(CC(l,d,i,j,k));

            for (int dir = 0; dir < 6; dir++)
            {
                const labelVector fo = faceOffsets[dir];

                if
                (
                    this->isInside(CC[l][d](ijk+fo))
                )
                {
                    wallAdjMask_(l,d,i,j,k) = 1;

                    // Neighbor cell in the IB
                    const vector nb(CC[l][d](ijk+fo));
                    const scalar wd = this->wallDistance(c,nb);
                    const scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                    wallDist_(l,d,i,j,k)[dir] = xi;
                }
            }
        }
    }
}

// Destructor

template<class Type, class MeshType>
Fadlun<Type,MeshType>::~Fadlun()
{}

template<class Type, class MeshType>
void Fadlun<Type,MeshType>::correctLinearSystem
(
    linearSystem<stencil,Type,MeshType>& ls
)
{
    forAllCells(ls.b(),l,d,i,j,k)
    {
        if (wallAdjMask_(l,d,i,j,k) == 1)
        {
            ls.b()(l,d,i,j,k) = Zero;
        }
    }

    forAllCells(ls.A(),l,d,i,j,k)
    {
        // Modify stencils in IB-adjacent cells
        if (wallAdjMask_(l,d,i,j,k) == 1)
        {
            scalar ximax = 0;
            // Loop over face number directions
            for (int dir = 0; dir < 6; dir++)
            {
                const label oppositeDir =
                    faceNumber(-faceOffsets[dir]);

                if (wallDist_(l,d,i,j,k)[dir] > ximax)
                {
                    ximax = wallDist_(l,d,i,j,k)[dir];
                    const scalar xic = 1.0 - wallDist_(l,d,i,j,k)[dir];
                    const scalar xinb = 1.0 + xic;
                    const scalar w = xic/xinb;

                    ls.A()(l,d,i,j,k) = Zero;
                    ls.A()(l,d,i,j,k).center() = 1.0;
                    ls.A()(l,d,i,j,k)[oppositeDir+1] = -w;
                }
            }
        }
    }
}

}

}

}
