#include "slipBoundaryCondition.H"

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
void slipBoundaryCondition<Type,colocated>::evaluate(const label d)
{
    const labelVector bo(this->offset());
    const label f(faceNumber(bo));
    const label l(this->l_);
    const label fd = f/2;

    meshDirection<Type,colocated>& field = this->level_[d];

    const colocatedVectorFaceField& fn = this->faceNormals();

    const labelVector S(this->S(d));
    const labelVector E(this->E(d));

    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        const labelVector upp(upperFaceNeighbor(ijk,f));

        field(ijk+bo) =
            field(ijk) - 2.0*(fn[fd](l,d,upp) & field(ijk))*fn[fd](l,d,upp);
    }
}

// Staggered

template<class Type>
void slipBoundaryCondition<Type,staggered>::eliminateGhosts
(
    linearSystem<stencil,Type,staggered>& sys,
    const label d
)
{
    const labelVector bo(this->offset());
    const label f(faceNumber(bo));
    const label l(this->l_);
    const label fd = f/2;

    if (fd == d)
    {
        meshField<stencil,staggered>& A = sys.A();
        meshField<Type,staggered>& b = sys.b();

        const meshField<scalar,staggered>& cv = this->cellVolumes();

        const labelVector S(this->S(d));
        const labelVector E(this->E(d));

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
        NeumannBoundaryCondition<Type,staggered>::eliminateGhosts(sys,d);
    }
}

template<class Type>
void slipBoundaryCondition<Type,staggered>::evaluate(const label d)
{
    const labelVector bo(this->offset());
    const label f(faceNumber(bo));
    const label fd = f/2;

    if (fd == d)
    {
        meshDirection<Type,staggered>& field = this->level_[d];

        const labelVector S(this->S(d));
        const labelVector E(this->E(d));

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
        NeumannBoundaryCondition<Type,staggered>::evaluate(d);
    }
}

}

}

}

