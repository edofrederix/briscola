#include "stagNeumannBoundaryCondition.H"

#include "colocated.H"
#include "staggered.H"
#include "meshLevel.H"
#include "linearSystem.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
stagNeumannBoundaryCondition<Type>::stagNeumannBoundaryCondition
(
    const meshField<Type,staggered>& mshField,
    const boundary& b,
    const List<Type>& boundaryGradients
)
:
    NeumannBoundaryCondition<Type,staggered>(mshField, b, boundaryGradients)
{}

template<class Type>
stagNeumannBoundaryCondition<Type>::stagNeumannBoundaryCondition
(
    const meshField<Type,staggered>& mshField,
    const boundary& b
)
:
    NeumannBoundaryCondition<Type,staggered>(mshField, b)
{}

template<class Type>
stagNeumannBoundaryCondition<Type>::stagNeumannBoundaryCondition
(
    const stagNeumannBoundaryCondition<Type>& bc
)
:
    NeumannBoundaryCondition<Type,staggered>(bc)
{}

template<class Type>
stagNeumannBoundaryCondition<Type>::stagNeumannBoundaryCondition
(
    const meshField<Type,staggered>& field,
    const stagNeumannBoundaryCondition<Type>& bc
)
:
    NeumannBoundaryCondition<Type,staggered>(field, bc)
{}

template<class Type>
void stagNeumannBoundaryCondition<Type>::eliminateGhosts
(
    linearSystem<stencil,Type,staggered>& sys,
    const label l,
    const label d
)
{
    const labelVector bo(this->offset());
    const label faceNum(faceNumber(bo));
    const label faceNum2(faceNumber(-bo));

    if (faceNum/2 == d)
    {
        meshField<stencil,staggered>& A = sys.A();
        meshField<Type,staggered>& b = sys.b();

        const meshField<faceScalar,staggered>& delta = this->faceDeltas();

        const labelVector S(this->S(l,d));
        const labelVector E(this->E(l,d));

        labelVector ijk;
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            const scalar ghostCoeff = A(l,d,ijk)[faceNum+1];

            A(l,d,ijk)[faceNum2+1] += ghostCoeff;
            A(l,d,ijk)[faceNum+1] = Zero;

            const scalar dx2 =
                1.0/delta(l,d,ijk)[faceNum] + 1.0/delta(l,d,ijk)[faceNum2];

            if (l == 0)
                b(l,d,ijk) -=
                    ghostCoeff*dx2*this->boundaryGradients_[d];
        }
    }
    else
    {
        NeumannBoundaryCondition<Type,staggered>::eliminateGhosts(sys,l,d);
    }
}

template<class Type>
void stagNeumannBoundaryCondition<Type>::eliminateGhosts
(
    linearSystem<symmStencil,Type,staggered>& sys,
    const label l,
    const label d
)
{
    NotImplemented;
}

template<class Type>
void stagNeumannBoundaryCondition<Type>::evaluate
(
    const label l,
    const label d
)
{
    const labelVector bo(this->offset());
    const label faceNum(faceNumber(bo));
    const label faceNum2(faceNumber(-bo));

    if (faceNum/2 == d)
    {
        const scalar H = l == 0;

        meshDirection<Type,staggered>& fd = this->mshField_[l][d];

        const meshDirection<faceScalar,staggered>& delta =
            this->faceDeltas()[l][d];

        const labelVector S(this->S(l,d));
        const labelVector E(this->E(l,d));

        labelVector ijk;
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            const scalar dx2 =
                1.0/delta(ijk)[faceNum] + 1.0/delta(ijk)[faceNum2];

            fd(ijk+bo) = fd(ijk-bo) + H*dx2*this->boundaryGradients_[d];
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

