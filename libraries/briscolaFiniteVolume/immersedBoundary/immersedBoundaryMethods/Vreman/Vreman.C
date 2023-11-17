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
    ),
    exchangePoints_(this->mask_.numberOfLevels())
{
    forAll(exchangePoints_, l)
    {
        exchangePoints_[l].setSize(MeshType::numberOfDirections);
    }

    // Cell centers
    const meshField<vector,MeshType>& CC =
        fvMsh.metrics<MeshType>().cellCenters();

    // Set IB mask fields
    forAllLevels(wallDist_,l,d,i,j,k)
    {
        const labelVector ijk(i,j,k);

        wallDist_(l,d,i,j,k) = -1.0;
        neighborDist_(l,d,i,j,k) = -1.0;

        if (this->ghostMask_(l,d,i,j,k) == 1)
        {
            // Ghost cell
            const vector gc(CC(l,d,i,j,k));

            for (int dir = 0; dir < 6; dir++)
            {
                const labelVector fo = faceOffsets[dir];

                if (!this->isInside(CC[l][d](ijk+fo)))
                {
                    // Wall-adjacent cell
                    const vector wa(CC[l][d](ijk+fo));
                    const scalar wd = this->wallDistance(wa,gc);
                    const scalar xi = (mag(gc-wa)-wd)/mag(gc-wa);

                    wallDist_(l,d,i,j,k)[dir] = xi;

                    // Second neighbor
                    const vector sn(CC[l][d](ijk+2.0*fo));
                    const scalar xi2 = mag(gc-sn)/mag(gc-wa);

                    neighborDist_(l,d,i,j,k)[dir] = xi2;

                    if
                    (
                           (i+2.0*fo.x() < 0 || i+2.0*fo.x() > neighborDist_[l][d].I().right())
                        || (j+2.0*fo.y() < 0 || j+2.0*fo.y() > neighborDist_[l][d].I().top())
                        || (k+2.0*fo.z() < 0 || k+2.0*fo.z() > neighborDist_[l][d].I().fore())
                    )
                    {
                        exchangePoints_[l][d].append(ijk+2.0*fo);
                    }

                    break;
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

    const fvMesh& fvMsh = x.fvMsh();

    label l = x.levelNum();

    forAll(x, d)
    {
        cellDataExchange<MeshType> exchange(exchangePoints_[l][d], fvMsh, l, d);

        List<Type> exchangeData(move(exchange(x.mshField())));

        scalar cursor = 0;

        forAllCells(x[d],i,j,k)
        {
            if (this->ghostMask_(l,d,i,j,k) == 1)
            {
                const labelVector ijk(i,j,k);
                for (int dir = 0; dir < 6; dir++)
                {
                    const labelVector fo = faceOffsets[dir];

                    if (this->wallDist_(l,d,i,j,k)[dir] > 0)
                    {
                        const scalar xi
                            = wallDist_(l,d,i,j,k)[dir];

                        scalar w1 = 2.0 - (2.0 - xi);
                        scalar w2 = -1.0 + (1.0 - xi);

                        Type secondNeighborValue;

                        if
                        (
                               (i+2.0*fo.x() < 0 || i+2.0*fo.x() > neighborDist_[l][d].I().right())
                            || (j+2.0*fo.y() < 0 || j+2.0*fo.y() > neighborDist_[l][d].I().top())
                            || (k+2.0*fo.z() < 0 || k+2.0*fo.z() > neighborDist_[l][d].I().fore())
                        )
                        {
                            secondNeighborValue = exchangeData[cursor++];
                        }
                        else
                        {
                            secondNeighborValue = x[d](ijk+2.0*fo);
                        }

                        x(d,i,j,k) = (1.0 - omega) * x(d,i,j,k)
                            + omega *
                            (
                                w1*x[d](ijk+fo) + w2*secondNeighborValue
                            );

                        break;
                    }
                }
            }
        }
    }


}

}

}

}
