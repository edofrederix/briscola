#include "mappedBoundaryCondition.H"
#include "mappedBoundary.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Colocated

template<class Type>
mappedBoundaryCondition<Type,colocated>::mappedBoundaryCondition
(
    const meshField<Type,colocated>& mshField,
    const boundary& b
)
:
    mappedBoundaryConditionBase<Type,colocated>(mshField,b)
{
    const vector mo(b.cast<mappedBoundary>().mappingOffset());
    const labelVector bo(this->offset());
    const label f = faceNumber(bo);

    // Set data exchanges

    this->exchanges_.clear();
    this->exchanges_.setSize(1);

    const meshDirection<faceVector,colocated>& fc =
        this->faceCenters()[0][0];

    const labelVector S(this->S(0,0));
    const labelVector E(this->E(0,0));
    const labelVector N(this->N(0,0));

    vectorList points(cmptProduct(N));

    // The sample points are taken as the boundary face centers offset by
    // the mapped boundary offset

    label c = 0;
    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        points[c++] = fc(ijk)[f] + mo;
    }

    this->exchanges_.set
    (
        0,
        new pointDataExchange<colocated>
        (
            points,
            this->fvMsh_,
            0,
            0
        )
    );
}

template<class Type>
void mappedBoundaryCondition<Type,colocated>::prepare(const label l)
{
    if (l > 0)
        return;

    const labelVector bo(this->offset());
    const label f = faceNumber(bo);

    meshField<Type,colocated>& field = this->mshField();

    pointDataExchange<colocated>& exchange = this->exchanges_[0];

    Field<Type> data(move(exchange(field)));

    const labelVector S(this->S(0,0));
    const labelVector E(this->E(0,0));

    if (this->setAverages_[0])
    {
        const meshDirection<faceScalar,colocated>& fa =
            this->faceAreas()[0][0];

        Type sum = Zero;
        scalar area = Zero;

        label c = 0;
        labelVector ijk;
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            const scalar a = fa(ijk-S)[f];

            area += a;
            sum += a*data[c++];
        }

        reduce(sum, sumOp<Type>());
        reduce(area, sumOp<scalar>());

        const Type average = sum/area;

        data += (this->averages_[0] - average);
    }

    block<Type>& bv = this->boundaryValues_[0];

    // Set the boundary values and let the Dirichlet boundary condition handle
    // the rest

    label c = 0;
    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        bv(ijk-S) = data[c++];
    }
}


// Staggered

template<class Type>
mappedBoundaryCondition<Type,staggered>::mappedBoundaryCondition
(
    const meshField<Type,staggered>& mshField,
    const boundary& b
)
:
    mappedBoundaryConditionBase<Type,staggered>(mshField,b)
{
    const vector mo(b.cast<mappedBoundary>().mappingOffset());
    const labelVector bo(this->offset());
    const label f = faceNumber(bo);

    // Set data exchanges

    this->exchanges_.clear();
    this->exchanges_.setSize(3);

    forAll(this->exchanges_, d)
    {
        const meshDirection<faceVector,staggered>& fc =
            this->faceCenters()[0][d];

        const meshDirection<vector,staggered>& cc =
            this->cellCenters()[0][d];

        const labelVector S(this->S(0,d));
        const labelVector E(this->E(0,d));
        const labelVector N(this->N(0,d));

        vectorList points(cmptProduct(N));

        // When the direction is shifted into the boundary then cell centers
        // coincide with the boundary and they are used as sample points.
        // Otherwise we take the face centers on the boundary as the sample
        // points.

        label c = 0;
        labelVector ijk;
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            points[c++] = (d == label(f/2) ? cc(ijk) : fc(ijk)[f]) + mo;
        }

        this->exchanges_.set
        (
            d,
            new pointDataExchange<staggered>
            (
                points,
                this->fvMsh_,
                0,
                d
            )
        );
    }
}

template<class Type>
void mappedBoundaryCondition<Type,staggered>::prepare(const label l)
{
    if (l > 0)
        return;

    const labelVector bo(this->offset());
    const label f = faceNumber(bo);

    meshField<Type,staggered>& field = this->mshField();

    forAll(field[l], d)
    {
        pointDataExchange<staggered>& exchange = this->exchanges_[d];

        Field<Type> data(move(exchange(field)));

        const labelVector S(this->S(0,d));
        const labelVector E(this->E(0,d));

        if (this->setAverages_[d])
        {
            // Colocated face areas needed for directions that are shifted into
            // the boundary
            const meshDirection<faceScalar,colocated>& cfa =
                this->fvMsh_.template metrics<colocated>().faceAreas()[0][0];

            // Staggered face areas
            const meshDirection<faceScalar,staggered>& sfa =
                this->faceAreas()[0][d];

            Type sum = Zero;
            scalar area = 0.0;

            label c = 0;
            labelVector ijk;
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                const scalar a = (d == label(f/2) ? cfa(ijk)[f] : sfa(ijk)[f]);

                area += a;
                sum += a*data[c++];
            }

            reduce(sum, sumOp<Type>());
            reduce(area, sumOp<scalar>());

            const Type average = sum/area;

            data += (this->averages_[d] - average);
        }

        block<Type>& bv = this->boundaryValues_[d];

        // Set the boundary values and let the Dirichlet boundary condition
        // handle the rest

        label c = 0;
        labelVector ijk;
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            bv(ijk-S) = data[c++];
        }
    }
}

}

}

}
