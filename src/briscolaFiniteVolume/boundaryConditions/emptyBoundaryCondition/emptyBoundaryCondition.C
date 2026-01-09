#include "emptyBoundaryCondition.H"

#include "colocated.H"
#include "staggered.H"
#include "meshFields.H"
#include "faceFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Colocated

template<class Type>
void emptyBoundaryCondition<Type,colocated>::eliminateGhosts
(
    linearSystem<stencil,Type,colocated>& sys,
    const label d
)
{
    const labelVector bo(this->offset());
    const label f(faceNumber(bo));
    const label l(this->l_);

    meshField<stencil,colocated>& A = sys.A();

    const labelVector S(this->S(d));
    const labelVector E(this->E(d));

    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        const scalar ghostCoeff = A(l,d,ijk)[f+1];

        A(l,d,ijk)[0] += ghostCoeff;
        A(l,d,ijk)[f+1] = Zero;
    }
}

template<class Type>
void emptyBoundaryCondition<Type,colocated>::evaluate(const label d)
{
    const labelVector bo(this->offset());

    meshDirection<Type,colocated>& field = this->level_[d];

    const labelVector S(this->S(d));
    const labelVector E(this->E(d));

    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        field(ijk+bo) = field(ijk);
    }
}

// Staggered

template<class Type>
void emptyBoundaryCondition<Type,staggered>::eliminateGhosts
(
    linearSystem<stencil,Type,staggered>& sys,
    const label d
)
{
    const labelVector bo(this->offset());
    const label f(faceNumber(bo));
    const label l(this->l_);
    const label fd = f/2;

    meshField<stencil,staggered>& A = sys.A();
    meshField<Type,staggered>& b = sys.b();

    const meshField<scalar,staggered>& cv = this->cellVolumes();

    const labelVector S(this->S(d));
    const labelVector E(this->E(d));

    labelVector ijk;

    if (fd == d)
    {

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
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            const scalar ghostCoeff = A(l,d,ijk)[f+1];

            A(l,d,ijk)[0] += ghostCoeff;
            A(l,d,ijk)[f+1] = Zero;
        }
    }
}

template<class Type>
void emptyBoundaryCondition<Type,staggered>::evaluate(const label d)
{
    const labelVector bo(this->offset());
    const label f(faceNumber(bo));
    const label fd = f/2;

    meshDirection<Type,staggered>& field = this->level_[d];

    const labelVector S(this->S(d));
    const labelVector E(this->E(d));

    labelVector ijk;

    if (fd == d)
    {
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            field(ijk+bo) = Zero;
        }
    }
    else
    {
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            field(ijk+bo) = field(ijk);
        }
    }
}

}

}

}

