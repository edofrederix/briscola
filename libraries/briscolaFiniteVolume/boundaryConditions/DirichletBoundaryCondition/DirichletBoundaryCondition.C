#include "DirichletBoundaryCondition.H"

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
DirichletBoundaryCondition<Type,MeshType>::DirichletBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b
)
:
    boundaryCondition<Type,MeshType>(mshField, b),
    boundaryValues_(this->dict().lookup("values"))
{}

template<class Type, class MeshType>
DirichletBoundaryCondition<Type,MeshType>::DirichletBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b,
    const List<Type>& boundaryValues
)
:
    boundaryCondition<Type,MeshType>(mshField, b),
    boundaryValues_(boundaryValues)
{}

template<class Type, class MeshType>
DirichletBoundaryCondition<Type,MeshType>::DirichletBoundaryCondition
(
    const DirichletBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(bc),
    boundaryValues_(bc.boundaryValues_)
{}

template<class Type, class MeshType>
DirichletBoundaryCondition<Type,MeshType>::DirichletBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const DirichletBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(field, bc),
    boundaryValues_(bc.boundaryValues_)
{}

template<class Type, class MeshType>
void DirichletBoundaryCondition<Type,MeshType>::eliminateGhosts
(
    linearSystem<stencil,Type,MeshType>& sys,
    const label l,
    const label d
)
{
    const labelVector bo(this->offset());
    const label faceNum(faceNumber(bo));

    meshField<stencil,MeshType>& A = sys.A();
    meshField<Type,MeshType>& b = sys.b();

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        const scalar ghostCoeff = A(l,d,ijk)[faceNum+1];

        A(l,d,ijk)[0] -= ghostCoeff;
        A(l,d,ijk)[faceNum+1] = Zero;

        if (l == 0)
            b(l,d,ijk) -= ghostCoeff*2.0*boundaryValues_[d];
    }
}

template<class Type, class MeshType>
void DirichletBoundaryCondition<Type,MeshType>::eliminateGhosts
(
    linearSystem<symmStencil,Type,MeshType>& sys,
    const label l,
    const label d
)
{
    const labelVector bo(this->offset());
    const label faceNum(faceNumber(bo));

    const labelVector o(faceNum%2 ? bo : zeroXYZ);

    meshField<symmStencil,MeshType>& A = sys.A();
    meshField<Type,MeshType>& b = sys.b();

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        const labelVector nei(ijk+o);
        const scalar ghostCoeff = A(l,d,nei)[faceNum/2+1];

        A(l,d,ijk)[0] -= ghostCoeff;
        A(l,d,nei)[faceNum/2+1] = Zero;

        if (l == 0)
            b(l,d,ijk) -= ghostCoeff*2.0*boundaryValues_[d];
    }
}

template<class Type, class MeshType>
void DirichletBoundaryCondition<Type,MeshType>::evaluate
(
    const label l,
    const label d
)
{
    const scalar H = l == 0;

    const labelVector bo(this->offset());

    meshDirection<Type,MeshType>& fd = this->mshField_[l][d];

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        fd(ijk+bo) = H*2.0*boundaryValues_[d] - fd(ijk);
    }
}

}

}

}

