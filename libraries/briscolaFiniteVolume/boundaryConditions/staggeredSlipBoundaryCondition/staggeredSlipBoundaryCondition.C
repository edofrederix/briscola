#include "staggeredSlipBoundaryCondition.H"

#include "colocated.H"
#include "staggered.H"

#include "meshLevel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
staggeredSlipBoundaryCondition<Type,MeshType>::staggeredSlipBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b
)
:
    boundaryCondition<Type,MeshType>(mshField, b)
{}

template<class Type, class MeshType>
staggeredSlipBoundaryCondition<Type,MeshType>::staggeredSlipBoundaryCondition
(
    const staggeredSlipBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(bc.mshField(), bc.mshBoundary())
{}

template<class Type, class MeshType>
staggeredSlipBoundaryCondition<Type,MeshType>::staggeredSlipBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const staggeredSlipBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(field, bc.mshBoundary())
{}

template<class Type, class MeshType>
void staggeredSlipBoundaryCondition<Type,MeshType>::prepare(const label)
{}

template<class Type, class MeshType>
void staggeredSlipBoundaryCondition<Type,MeshType>::evaluate(const label l)
{
    meshLevel<Type,MeshType>& field = this->mshField()[l];

    const labelVector bo(this->offset());

    forAll(field, d)
    {
        meshDirection<Type,MeshType>& fd = field[d];

        const labelVector S(this->S(l,d));
        const labelVector E(this->E(l,d));

        labelVector ijk;

        // Zero value on shifted boundaries, zero gradient otherwise

        if (MeshType::shifted(d,bo))
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                fd(ijk+bo) = Zero;
            }
        }
        else
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                fd(ijk+bo) = fd(ijk);
            }
        }
    }
}

template<class Type, class MeshType>
tmp<block<Type>> staggeredSlipBoundaryCondition<Type,MeshType>::internalValue
(
    const label l,
    const label d
)
{
    return tmp<block<Type>>
    (
        new block<Type>(this->N(l,d), Zero)
    );
}

template<class Type, class MeshType>
tmp<block<Type>> staggeredSlipBoundaryCondition<Type,MeshType>::boundarySources
(
    const label l,
    const label d
)
{
    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    tmp<block<Type>> tR(new block<Type>(E-S), Zero);

    return tR;
}

}

}

}

