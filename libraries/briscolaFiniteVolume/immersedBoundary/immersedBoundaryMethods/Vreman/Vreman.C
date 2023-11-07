#include "Vreman.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

template<class Type, class MeshType>
Vreman<Type,MeshType>::Vreman
(
    dictionary& dict,
    const fvMesh& fvMsh
)
:
    immersedBoundaryMethod<Type,MeshType>(dict,fvMsh,true),
    mask_
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
    ),
    neighborDist_
    (
        "neighborDist",
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
    forAllLevels(mask_,l,d,i,j,k)
    {
        const labelVector ijk(i,j,k);

        mask_(l,d,i,j,k) = 0;
        wallDist_(l,d,i,j,k) = -1.0;
        neighborDist_(l,d,i,j,k) = -1.0;

        if (this->isInside(CC(l,d,i,j,k)))
        {
            mask_(l,d,i,j,k) = 1;
        }
        else
        {
            const vector c(CC(l,d,i,j,k));

            for (int dir = 0; dir < 6; dir++)
            {
                const labelVector fo = faceOffsets[dir];
                const label oppositeDir = faceNumber(-fo);

                if
                (
                    this->isInside(CC[l][d](ijk+fo))
                )
                {
                    // Neighbor cell in the IB
                    const vector nb(CC[l][d](ijk+fo));
                    const scalar wd = this->wallDistance(c,nb);
                    const scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                    wallDist_(l,d,i,j,k)[dir] = xi;

                    // Neighbor cell in the opposite direction
                    const vector nbo(CC[l][d](ijk-fo));
                    const scalar xi2 = mag(nbo-nb)/mag(c-nb);

                    neighborDist_(l,d,i,j,k)[oppositeDir] = xi2;
                }
            }
        }
    }
}

// Destructor

template<class Type, class MeshType>
Vreman<Type,MeshType>::~Vreman()
{}

template<class Type, class MeshType>
void Vreman<Type,MeshType>::correctJacobiPoints
(
    meshLevel<Type,MeshType>& x
)
{
    scalar omega = this->omega_;
    label l = x.levelNum();

    forAllDirections(x,d,i,j,k)
    {
        const labelVector ijk(i,j,k);
        if (this->ghostMask_(l,d,i,j,k) == 1)
        {
            for (int dir = 0; dir < 6; dir++)
            {
                const label oppositeDir =
                    faceNumber(-faceOffsets[dir]);

                const labelVector neighbor = ijk + faceOffsets[dir];
                const labelVector secondNeighbor
                    = ijk + 2.0*faceOffsets[dir];

                if (mask_[l][d](ijk+faceOffsets[dir]) == 0)
                {
                    const scalar xi
                        = wallDist_[l][d](neighbor)[oppositeDir];

                    scalar w1 = 2.0 - (2.0 - xi);
                    scalar w2 = -1.0 + (1.0 - xi);

                    // if second neighbor is inside processor
                    x(d,i,j,k) = (1.0 - omega) * x(d,i,j,k)
                        + omega *
                        (
                            w1*x[d](neighbor) + w2*x[d](secondNeighbor)
                        );

                    // else { ... }
                    break;
                }
            }
        }
    }
}

}

}

}
