#include "NeumannBoundaryCondition.H"

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
NeumannBoundaryCondition<Type,MeshType>::NeumannBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b
)
:
    boundaryCondition<Type,MeshType>(mshField, b),
    boundaryGradients_(this->dict().lookup("gradients"))
{}

template<class Type, class MeshType>
NeumannBoundaryCondition<Type,MeshType>::NeumannBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b,
    const List<Type>& boundaryGradients
)
:
    boundaryCondition<Type,MeshType>(mshField, b),
    boundaryGradients_(boundaryGradients)
{}

template<class Type, class MeshType>
NeumannBoundaryCondition<Type,MeshType>::NeumannBoundaryCondition
(
    const NeumannBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(bc),
    boundaryGradients_(bc.boundaryGradients_)
{}

template<class Type, class MeshType>
NeumannBoundaryCondition<Type,MeshType>::NeumannBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const NeumannBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(field, bc),
    boundaryGradients_(bc.boundaryGradients_)
{}

template<class Type, class MeshType>
void NeumannBoundaryCondition<Type,MeshType>::eliminateGhosts
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

    const meshField<faceScalar,MeshType>& delta = this->faceDeltas();

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        const scalar ghostCoeff = A(l,d,ijk)[faceNum+1];

        A(l,d,ijk)[0] += ghostCoeff;
        A(l,d,ijk)[faceNum+1] = Zero;

        if (l == 0)
        {
            const scalar dx = 1.0/delta(l,d,ijk)[faceNum];
            b(l,d,ijk) -= ghostCoeff*dx*boundaryGradients_[d];
        }
    }
}

template<class Type, class MeshType>
void NeumannBoundaryCondition<Type,MeshType>::evaluate
(
    const label l,
    const label d
)
{
    const scalar H = l == 0;

    const labelVector bo(this->offset());
    const label faceNum(faceNumber(bo));

    meshDirection<Type,MeshType>& fd = this->mshField_[l][d];

    const meshDirection<faceScalar,MeshType>& delta =
        this->faceDeltas()[l][d];

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        const scalar dx = 1.0/delta(ijk)[faceNum];
        fd(ijk+bo) = fd(ijk) + H*dx*boundaryGradients_[d];
    }
}

}

}

}

