#include "emptyBoundaryCondition.H"

#include "colocated.H"
#include "staggered.H"
#include "meshField.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Staggered

template<class Type>
void emptyBoundaryCondition<Type,staggered>::eliminateGhosts
(
    linearSystem<stencil,Type,staggered>& sys,
    const label l,
    const label d
)
{
    const labelVector bo(this->offset());
    const label f(faceNumber(bo));
    const label fd = f/2;

    if (fd == d)
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
        emptyBoundaryConditionBase<Type,staggered>::eliminateGhosts(sys,l,d);
    }
}

template<class Type>
void emptyBoundaryCondition<Type,staggered>::evaluate
(
    const label l,
    const label d
)
{
    const labelVector bo(this->offset());
    const label f(faceNumber(bo));
    const label fd = f/2;

    if (fd == d)
    {
        meshDirection<Type,staggered>& field = this->mshField_[l][d];

        const labelVector S(this->S(l,d));
        const labelVector E(this->E(l,d));

        labelVector ijk;
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            field(ijk+bo) = Zero;
        }
    }
    else
    {
        emptyBoundaryConditionBase<Type,staggered>::evaluate(l,d);
    }
}

}

}

}

