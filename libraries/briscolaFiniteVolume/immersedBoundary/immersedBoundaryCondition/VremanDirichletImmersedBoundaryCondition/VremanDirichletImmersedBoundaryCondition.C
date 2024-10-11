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
    exchangePoints_(MeshType::numberOfDirections),
    boundaryValues_(this->dict().lookup("values"))
{
    // Check shape overlap
    if (this->IB_.shapeOverlap())
    {
        WarningInFunction
            << "Overlapping shapes identified."
            << " This may cause issues with Vreman IBM." << endl;
    }

    // Cell centers
    const meshField<vector,MeshType>& CC =
        this->fvMsh_.template metrics<MeshType>().cellCenters();

    // Set exchange points list
    forAllCells(this->forcingPoints_[0],d,i,j,k)
    {
        const labelVector ijk(i,j,k);

        if (this->forcingPoints_(0,d,i,j,k))
        {
            faceLabel I = this->fvMsh_.template I<MeshType>(0,d);

            for (int dir = 0; dir < 6; dir++)
            {
                const labelVector fo = faceOffsets[dir];

                if (!this->IB_.mask()[0][d](ijk+fo))
                {
                    if
                    (
                           (i+2.0*fo.x() < I.left())
                        || (i+2.0*fo.x() > I.right())
                        || (j+2.0*fo.y() < I.bottom())
                        || (j+2.0*fo.y() > I.top())
                        || (k+2.0*fo.z() < I.aft())
                        || (k+2.0*fo.z() > I.fore())
                    )
                    {
                        exchangePoints_[d].append(ijk+2.0*fo);
                    }

                    if (this->IB_.isInside(CC[0][d](ijk+2.0*fo)))
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

    if (l == 0)
    {
        forAll(x, d)
        {
            faceLabel I = fvMsh.template I<MeshType>(0,d);

            cellDataExchange<MeshType> exchange
            (
                exchangePoints_[d],
                fvMsh,
                l,
                d
            );

            List<Type> exchangeData(move(exchange(x.mshField())));

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
                                   (i+2.0*fo.x() < I.left())
                                || (i+2.0*fo.x() > I.right())
                                || (j+2.0*fo.y() < I.bottom())
                                || (j+2.0*fo.y() > I.top())
                                || (k+2.0*fo.z() < I.aft())
                                || (k+2.0*fo.z() > I.fore())
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

}
