#include "VremanDirichletImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

template<class Type, class MeshType>
VremanDirichletImmersedBoundaryCondition<Type,MeshType>::
VremanDirichletImmersedBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const immersedBoundary<MeshType>& ib
)
:
    immersedBoundaryCondition<Type,MeshType>
    (
        mshField,
        ib,
        ib.ghostMask()
    ),
    exchangePoints_(this->IB_.mask().numberOfLevels()),
    exchanges_(this->IB_.mask().numberOfLevels()),
    boundaryValues_(this->dict().lookup("values"))
{
    // Check shape overlap
    if (this->IB_.shapeOverlap())
    {
        WarningInFunction
            << "Overlapping shapes identified."
            << " This may cause issues with Vreman IBM." << endl;
    }

    forAll(exchangePoints_, l)
    {
        exchangePoints_[l].setSize(MeshType::numberOfDirections);
        exchanges_[l].setSize(MeshType::numberOfDirections);
    }

    // Cell centers
    const meshField<vector,MeshType>& CC =
        this->fvMsh_.template metrics<MeshType>().cellCenters();

    // Set exchange points list
    forAllCells(this->forcingPoints_,l,d,i,j,k)
    {
        const labelVector ijk(i,j,k);

        if (this->forcingPoints_(l,d,i,j,k))
        {
            faceLabel I = this->fvMsh_.template I<MeshType>(l,d);

            for (int dir = 0; dir < 6; dir++)
            {
                const labelVector fo = faceOffsets[dir];

                if (!this->IB_.mask()[l][d](ijk+fo))
                {
                    if
                    (
                           (i+2*fo.x() < I.left())
                        || (i+2*fo.x() > I.right())
                        || (j+2*fo.y() < I.bottom())
                        || (j+2*fo.y() > I.top())
                        || (k+2*fo.z() < I.aft())
                        || (k+2*fo.z() > I.fore())
                    )
                    {
                        exchangePoints_[l][d].append(ijk+2*fo);
                    }

                    if (this->IB_.isInside(CC[l][d](ijk+2*fo)))
                    {
                        WarningInFunction
                            << "Second neighbor point of ghost cell "
                            << vector(i,j,k) << " at d = "
                            << d << " located inside immersed boundary."
                            << " This may cause issues with Vreman IBM."
                            << endl;
                    }

                    break;
                }
            }
        }
    }

    // Set cell data exchanges
    forAll(exchanges_, l)
    {
        forAll(exchanges_[l], d)
        {
            exchanges_[l].set
            (
                d,
                new cellDataExchange<MeshType>
                (
                    exchangePoints_[l][d],
                    this->fvMsh_,
                    l,
                    d
                )
            );
        }
    }
}

// Destructor

template<class Type, class MeshType>
VremanDirichletImmersedBoundaryCondition<Type,MeshType>::
~VremanDirichletImmersedBoundaryCondition()
{}

template<class Type, class MeshType>
void VremanDirichletImmersedBoundaryCondition<Type,MeshType>::
correctJacobiPoints
(
    meshLevel<Type,MeshType>& x
) const
{
    scalar omega = this->omega_;

    const fvMesh& fvMsh = x.fvMsh();

    label l = x.levelNum();

    forAll(x, d)
    {
        faceLabel I = fvMsh.template I<MeshType>(l,d);

        List<Type> exchangeData
        (
            move(exchanges_[l][d].dataFunc(x.mshField()))
        );

        scalar cursor = 0;

        forAllCells(x[d],i,j,k)
        {
            if (this->forcingPoints_(l,d,i,j,k))
            {
                const labelVector ijk(i,j,k);

                for (int dir = 0; dir < 6; dir++)
                {
                    const labelVector fo = faceOffsets[dir];

                    if (this->IB_.wallDistGhost()(l,d,i,j,k)[dir] > 0)
                    {
                        const scalar xi
                            = this->IB_.wallDistGhost()(l,d,i,j,k)[dir];

                        scalar w1 = 2.0 - (2.0 - xi);
                        scalar w2 = -1.0 + (1.0 - xi);

                        Type secondNeighborValue;

                        if
                        (
                               (i+2*fo.x() < I.left())
                            || (i+2*fo.x() > I.right())
                            || (j+2*fo.y() < I.bottom())
                            || (j+2*fo.y() > I.top())
                            || (k+2*fo.z() < I.aft())
                            || (k+2*fo.z() > I.fore())
                        )
                        {
                            secondNeighborValue = exchangeData[cursor++];
                        }
                        else
                        {
                            secondNeighborValue = x[d](ijk+2*fo);
                        }

                        x(d,i,j,k) = (1.0 - omega) * x(d,i,j,k)
                            + omega *
                            (
                                boundaryValues_[d]
                                + w1*x[d](ijk+fo)
                                + w2*secondNeighborValue
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
