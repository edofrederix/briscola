#include "stagDirichletBoundaryCondition.H"

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
stagDirichletBoundaryCondition<Type>::stagDirichletBoundaryCondition
(
    const meshField<Type,staggered>& mshField,
    const boundary& b
)
:
    DirichletBoundaryCondition<Type,staggered>(mshField, b)
{}

template<class Type>
stagDirichletBoundaryCondition<Type>::stagDirichletBoundaryCondition
(
    const meshField<Type,staggered>& mshField,
    const boundary& b,
    const List<Type>& boundaryValues
)
:
    DirichletBoundaryCondition<Type,staggered>(mshField, b, boundaryValues)
{}

template<class Type>
stagDirichletBoundaryCondition<Type>::stagDirichletBoundaryCondition
(
    const stagDirichletBoundaryCondition<Type>& bc
)
:
    DirichletBoundaryCondition<Type,staggered>(bc)
{}

template<class Type>
stagDirichletBoundaryCondition<Type>::stagDirichletBoundaryCondition
(
    const meshField<Type,staggered>& field,
    const stagDirichletBoundaryCondition<Type>& bc
)
:
    DirichletBoundaryCondition<Type,staggered>(field, bc)
{}

template<class Type>
void stagDirichletBoundaryCondition<Type>::eliminateGhosts
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
                b(l,d,ijk) = this->boundaryValues_[d]*cv(l,d,ijk);
        }
    }
    else
    {
        DirichletBoundaryCondition<Type,staggered>::eliminateGhosts(sys,l,d);
    }
}

template<class Type>
void stagDirichletBoundaryCondition<Type>::evaluate
(
    const label l,
    const label d
)
{

    const labelVector bo(this->offset());
    const label faceNum(faceNumber(bo));

    if (faceNum/2 == d)
    {
        const scalar H = l == 0;

        meshDirection<Type,staggered>& fd = this->mshField_[l][d];

        const labelVector S(this->S(l,d));
        const labelVector E(this->E(l,d));

        labelVector ijk;
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            fd(ijk+bo) = H*this->boundaryValues_[d];
        }
    }
    else
    {
        DirichletBoundaryCondition<Type,staggered>::evaluate(l,d);
    }
}

}

}

}

