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

// Colocated

template<class Type>
void DirichletBoundaryCondition<Type,colocated>::eliminateGhosts
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
            b(l,d,ijk) -= ghostCoeff*2.0*this->boundaryValues_[l](ijk-S);
    }
}

template<class Type>
void DirichletBoundaryCondition<Type,colocated>::evaluate
(
    const label l,
    const label d
)
{
    const scalar H = l == 0;

    const labelVector bo(this->offset());

    meshDirection<Type,colocated>& field = this->mshField_[l][d];

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        field(ijk+bo) =
            H*2.0*this->boundaryValues_[l](ijk-S) - field(ijk);
    }
}

// Staggered

template<class Type>
void DirichletBoundaryCondition<Type,staggered>::eliminateGhosts
(
    linearSystem<stencil,Type,staggered>& sys,
    const label l,
    const label d
)
{
    const label c = l*3 + d;

    const labelVector bo(this->offset());
    const label f(faceNumber(bo));
    const label fd = f/2;

    meshField<stencil,staggered>& A = sys.A();
    meshField<Type,staggered>& b = sys.b();

    const meshField<scalar,staggered>& cv = this->cellVolumes();

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;

    if (fd == d)
    {
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            A(l,d,ijk) = diagStencil(cv(l,d,ijk));

            if (l == 0)
                b(l,d,ijk) = this->boundaryValues_[c](ijk-S)*cv(l,d,ijk);
        }
    }
    else
    {
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            const scalar ghostCoeff = A(l,d,ijk)[f+1];

            A(l,d,ijk)[0] -= ghostCoeff;
            A(l,d,ijk)[f+1] = Zero;

            if (l == 0)
                b(l,d,ijk) -= ghostCoeff*2.0*this->boundaryValues_[c](ijk-S);
        }
    }
}

template<class Type>
void DirichletBoundaryCondition<Type,staggered>::evaluate
(
    const label l,
    const label d
)
{
    const label c = l*3 + d;

    const labelVector bo(this->offset());
    const label f(faceNumber(bo));
    const label fd = f/2;

    const scalar H = l == 0;

    meshDirection<Type,staggered>& field = this->mshField_[l][d];

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;

    if (fd == d)
    {
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            field(ijk+bo) = H*this->boundaryValues_[c](ijk-S);
        }
    }
    else
    {
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            field(ijk+bo) =
                H*2.0*this->boundaryValues_[c](ijk-S) - field(ijk);
        }
    }
}

}

}

}

