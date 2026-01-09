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
    const meshLevel<Type,colocated>& level,
    const boundary& b
)
:
    mappedBoundaryConditionBase<Type,colocated>(level,b)
{
    const vector mo(b.cast<mappedBoundary>().mappingOffset());
    const labelVector bo(this->offset());
    const label f = faceNumber(bo);
    const label fd = f/2;

    // Set data exchanges

    this->exchanges_.clear();
    this->exchanges_.setSize(1);

    const faceField<vector,colocated>& fc = this->faceCenters();

    const labelVector S(this->S(0));
    const labelVector E(this->E(0));
    const labelVector N(this->N(0));

    vectorList points(cmptProduct(N));

    // The sample points are taken as the boundary face centers offset by
    // the mapped boundary offset

    label c = 0;
    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        points[c++] = fc[fd](upperFaceNeighbor(ijk,f)) + mo;
    }

    this->exchanges_.set
    (
        0,
        new pointDataExchange<colocated>
        (
            points,
            this->fvMsh_,
            this->l_,
            0
        )
    );
}

template<class Type>
void mappedBoundaryCondition<Type,colocated>::prepare()
{
    const labelVector bo(this->offset());
    const label f = faceNumber(bo);
    const label l = this->l_;
    const label fd = f/2;

    if (l > 0)
        return;

    meshField<Type,colocated>& field = this->level_.field();

    pointDataExchange<colocated>& exchange = this->exchanges_[0];

    Field<Type> data(move(exchange(field)));

    const labelVector S(this->S(0));
    const labelVector E(this->E(0));

    if (this->setAverage_)
    {
        const faceField<scalar,colocated>& fa = this->faceAreas();

        Type sum = Zero;
        scalar area = Zero;

        label c = 0;
        labelVector ijk;
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            const scalar a = fa[fd](0,0,upperFaceNeighbor(ijk,f));

            area += a;
            sum += a*data[c++];
        }

        reduce(sum, sumOp<Type>());
        reduce(area, sumOp<scalar>());

        const Type average = sum/area;

        data += (this->average_ - average);
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
    const meshLevel<Type,staggered>& level,
    const boundary& b
)
:
    mappedBoundaryConditionBase<Type,staggered>(level,b)
{
    const vector mo(b.cast<mappedBoundary>().mappingOffset());
    const labelVector bo(this->offset());
    const label f = faceNumber(bo);
    const label fd = f/2;

    // Set data exchanges

    this->exchanges_.clear();
    this->exchanges_.setSize(3);

    forAll(this->exchanges_, d)
    {
        const faceField<vector,staggered>& fc = this->faceCenters();

        const meshDirection<vector,staggered>& cc =
            this->cellCenters()[0][d];

        const labelVector S(this->S(d));
        const labelVector E(this->E(d));
        const labelVector N(this->N(d));

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
            points[c++] =
                (d == fd ? cc(ijk) : fc[fd](0,d,upperFaceNeighbor(ijk,f))) + mo;
        }

        this->exchanges_.set
        (
            d,
            new pointDataExchange<staggered>
            (
                points,
                this->fvMsh_,
                this->l_,
                d
            )
        );
    }
}

template<class Type>
void mappedBoundaryCondition<Type,staggered>::prepare()
{
    const labelVector bo(this->offset());
    const label f = faceNumber(bo);
    const label l = this->l_;
    const label fd = f/2;

    if (l > 0)
        return;

    meshField<Type,staggered>& field = this->level_.field();

    const tensor base =
        this->fvMsh_.msh().template cast<rectilinearMesh>().base();

    forAll(field[l], d)
    {
        pointDataExchange<staggered>& exchange = this->exchanges_[d];

        Field<Type> data(move(exchange(field)));

        const labelVector S(this->S(d));
        const labelVector E(this->E(d));

        if (this->setAverage_)
        {
            // Colocated face areas needed for directions that are shifted into
            // the boundary
            const colocatedScalarFaceField& cfa =
                this->fvMsh_.template metrics<colocated>().faceAreas();

            // Staggered face areas
            const staggeredScalarFaceField& sfa = this->faceAreas();

            Type sum = Zero;
            scalar area = 0.0;

            label c = 0;
            labelVector ijk;
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                const labelVector upp(upperFaceNeighbor(ijk,f));

                const scalar a =
                    (d == fd ? cfa[fd](0,0,upp) : sfa[fd](0,d,upp));

                area += a;
                sum += a*data[c++];
            }

            reduce(sum, sumOp<Type>());
            reduce(area, sumOp<scalar>());

            const Type average = sum/area;

            data += (staggered::project(this->average_, d, base) - average);
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
