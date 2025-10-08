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

// Colocated

template<class Type>
void NeumannBoundaryCondition<Type,colocated>::eliminateGhosts
(
    linearSystem<stencil,Type,colocated>& sys,
    const label l,
    const label d
)
{
    const labelVector bo(this->offset());
    const label faceNum(faceNumber(bo));

    meshField<stencil,colocated>& A = sys.A();
    meshField<Type,colocated>& b = sys.b();

    const meshField<faceScalar,colocated>& delta = this->faceDeltas();

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
            b(l,d,ijk) -= ghostCoeff*dx*this->boundaryGradients_[l](ijk-S);
        }
    }
}

template<class Type>
void NeumannBoundaryCondition<Type,colocated>::evaluate
(
    const label l,
    const label d
)
{
    const scalar H = l == 0;

    const labelVector bo(this->offset());
    const label faceNum(faceNumber(bo));

    meshDirection<Type,colocated>& fd = this->mshField_[l][d];

    const meshDirection<faceScalar,colocated>& delta =
        this->faceDeltas()[l][d];

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        const scalar dx = 1.0/delta(ijk)[faceNum];
        fd(ijk+bo) = fd(ijk) + H*dx*this->boundaryGradients_[l](ijk-S);
    }
}

// Staggered

template<class Type>
void NeumannBoundaryCondition<Type,staggered>::eliminateGhosts
(
    linearSystem<stencil,Type,staggered>& sys,
    const label l,
    const label d
)
{
    const label c = l*3 + d;

    const labelVector bo(this->offset());
    const label faceNum(faceNumber(bo));
    const label faceNum2(faceNumber(-bo));

    meshField<stencil,staggered>& A = sys.A();
    meshField<Type,staggered>& b = sys.b();

    const meshField<faceScalar,staggered>& delta = this->faceDeltas();

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;

    if (faceNum/2 == d)
    {
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
                    ghostCoeff*dx2*this->boundaryGradients_[c](ijk-S);
        }
    }
    else
    {
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
                b(l,d,ijk) -= ghostCoeff*dx*this->boundaryGradients_[c](ijk-S);
            }
        }
    }
}

template<class Type>
void NeumannBoundaryCondition<Type,staggered>::evaluate
(
    const label l,
    const label d
)
{
    const label c = l*3 + d;

    const labelVector bo(this->offset());
    const label faceNum(faceNumber(bo));
    const label faceNum2(faceNumber(-bo));

    const scalar H = l == 0;

    meshDirection<Type,staggered>& fd = this->mshField_[l][d];

    const meshDirection<faceScalar,staggered>& delta =
        this->faceDeltas()[l][d];

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;

    if (faceNum/2 == d)
    {
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            const scalar dx2 =
                1.0/delta(ijk)[faceNum] + 1.0/delta(ijk)[faceNum2];

            fd(ijk+bo) = fd(ijk-bo) + H*dx2*this->boundaryGradients_[c](ijk-S);
        }
    }
    else
    {
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            const scalar dx = 1.0/delta(ijk)[faceNum];
            fd(ijk+bo) = fd(ijk) + H*dx*this->boundaryGradients_[c](ijk-S);
        }
    }
}

}

}

}

