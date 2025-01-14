#include "stagEmptyBoundaryCondition.H"

#include "colocated.H"
#include "staggered.H"
#include "meshLevel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
stagEmptyBoundaryCondition<Type>::stagEmptyBoundaryCondition
(
    const meshField<Type,staggered>& mshField,
    const boundary& b
)
:
    emptyBoundaryCondition<Type,staggered>(mshField, b)
{}

template<class Type>
stagEmptyBoundaryCondition<Type>::stagEmptyBoundaryCondition
(
    const stagEmptyBoundaryCondition<Type>& bc
)
:
    emptyBoundaryCondition<Type,staggered>(bc)
{}

template<class Type>
stagEmptyBoundaryCondition<Type>::stagEmptyBoundaryCondition
(
    const meshField<Type,staggered>& field,
    const stagEmptyBoundaryCondition<Type>& bc
)
:
    emptyBoundaryCondition<Type,staggered>(field, bc)
{}

template<class Type>
void stagEmptyBoundaryCondition<Type>::eliminateGhosts
(
    linearSystem<stencil,Type,staggered>& sys,
    const label l,
    const label d
)
{
    const labelVector bo(this->offset());
    const label faceNum(faceNumber(bo));

    if (faceNum/2 == d)
    {
        meshField<stencil,staggered>& A = sys.A();
        meshField<Type,staggered>& b = sys.b();

        const meshField<scalar,staggered>& cv = this->cellVolumes();

        const labelVector S(this->S(l,d));
        const labelVector E(this->E(l,d));

        labelVector ijk;
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            A(l,d,ijk) = diagStencil(cv(l,d,ijk));

            if (l == 0)
                b(l,d,ijk) = Zero;
        }
    }
    else
    {
        NeumannBoundaryCondition<Type,staggered>::eliminateGhosts(sys,l,d);
    }
}

template<class Type>
void stagEmptyBoundaryCondition<Type>::evaluate
(
    const label l,
    const label d
)
{

    const labelVector bo(this->offset());
    const label faceNum(faceNumber(bo));

    if (faceNum/2 == d)
    {
        meshDirection<Type,staggered>& fd = this->mshField_[l][d];

        const labelVector S(this->S(l,d));
        const labelVector E(this->E(l,d));

        labelVector ijk;
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            fd(ijk+bo) = Zero;
        }
    }
    else
    {
        NeumannBoundaryCondition<Type,staggered>::evaluate(l,d);
    }
}

}

}

}

