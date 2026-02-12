#include "VremanDirichletImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructors

template<class Type, class MeshType>
VremanDirichletImmersedBoundaryCondition<Type,MeshType>::
VremanDirichletImmersedBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const immersedBoundary<MeshType>& ib
)
:
    immersedBoundaryCondition<Type,MeshType>(field, ib, &ib.ghostMask()),
    exchanges_(field.msh().size()),
    boundaryValues_(this->read("value"))
{
    // Check shape overlap

    if (this->ib_.shapeOverlap())
        this->shapeOverlapWarning();

    // Set mirror points

    const meshField<label,MeshType>& mask = this->forcingMask();

    forAll(mask, l)
    {
        exchanges_.set
        (
            l,
            new PtrList<cellDataExchange<MeshType>>
            (
                MeshType::numberOfDirections
            )
        );

        forAll(mask[l], d)
        {
            label nPoints = 0;

            forAllCells(mask[l][d],i,j,k)
                if (mask(l,d,i,j,k))
                    nPoints++;

            List<labelVector> exchangePoints(nPoints);

            label c = 0;

            forAllCells(mask[l][d],i,j,k)
            if (mask(l,d,i,j,k))
            {
                const labelVector ijk(i,j,k);

                exchangePoints[c] = ijk;

                for (int f = 0; f < 6; f++)
                {
                    const labelVector fo = faceOffsets[f];

                    if (!this->ib_.mask()[l][d](ijk + fo))
                    {
                        exchangePoints[c] += 2*fo;
                        break;
                    }
                }

                if (exchangePoints[c] == ijk)
                    WarningInFunction
                        << "No valid neighbor point found." << endl;

                c++;
            }

            exchanges_[l].set
            (
                d,
                new cellDataExchange<MeshType>
                (
                    exchangePoints,
                    this->fvMsh_,
                    l,
                    d
                )
            );
        }
    }
}

template<class Type, class MeshType>
VremanDirichletImmersedBoundaryCondition<Type,MeshType>::
VremanDirichletImmersedBoundaryCondition
(
    const VremanDirichletImmersedBoundaryCondition<Type,MeshType>& ibc
)
:
    immersedBoundaryCondition<Type,MeshType>(ibc),
    exchanges_(ibc.exchanges_.size()),
    boundaryValues_(ibc.boundaryValues_)
{
    forAll(exchanges_, i)
    {
        PtrList<cellDataExchange<MeshType>> copy(ibc.exchanges_[i]);
        exchanges_[i].transfer(copy);
    }
}

template<class Type, class MeshType>
VremanDirichletImmersedBoundaryCondition<Type,MeshType>::
VremanDirichletImmersedBoundaryCondition
(
    const VremanDirichletImmersedBoundaryCondition<Type,MeshType>& ibc,
    const meshField<Type,MeshType>& field
)
:
    immersedBoundaryCondition<Type,MeshType>(ibc, field),
    exchanges_(ibc.exchanges_.size()),
    boundaryValues_(ibc.boundaryValues_)
{
    forAll(exchanges_, i)
    {
        PtrList<cellDataExchange<MeshType>> copy(ibc.exchanges_[i]);
        exchanges_[i].transfer(copy);
    }
}

// Destructor

template<class Type, class MeshType>
VremanDirichletImmersedBoundaryCondition<Type,MeshType>::
~VremanDirichletImmersedBoundaryCondition()
{}

template<class Type, class MeshType>
void VremanDirichletImmersedBoundaryCondition<Type,MeshType>::evaluate
(
    const label l,
    const label d
)
{
    const scalar omega = this->omega_;

    meshDirection<Type,MeshType>& x = this->field_[l][d];
    const meshDirection<label,MeshType>& mask = this->forcingMask()[l][d];
    const meshDirection<faceScalar,MeshType>& y =
        this->ib_.wallDistGhost()[l][d];

    List<Type> data(move(exchanges_[l][d](this->field_)));

    label c = 0;

    forAllCells(x,i,j,k)
    {
        const labelVector ijk(i,j,k);

        if (mask(ijk))
        {
            for (int f = 0; f < 6; f++)
            {
                const labelVector fo = faceOffsets[f];

                if (!this->ib_.mask()(ijk + fo))
                {
                    const scalar xi = y(ijk)[f];

                    scalar w1 =  2.0 - (2.0 - xi);
                    scalar w2 = -1.0 + (1.0 - xi);

                    x(ijk) =
                        (1.0 - omega)*x(ijk)
                      + omega
                      * (
                            boundaryValues_[d]
                          + w1*x(ijk + fo)
                          + w2*data[c]
                        );

                    break;
                }
            }

            c++;
        }
    }
}

}

}

}
