#include "stagSlipBoundaryCondition.H"

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
stagSlipBoundaryCondition<Type>::stagSlipBoundaryCondition
(
    const meshField<Type,staggered>& mshField,
    const boundary& b
)
:
    stagNeumannBoundaryCondition<Type>
    (
        mshField,
        b,
        List<Type>(staggered::numberOfDirections, Zero)
    )
{}

template<class Type>
stagSlipBoundaryCondition<Type>::stagSlipBoundaryCondition
(
    const stagSlipBoundaryCondition<Type>& bc
)
:
    stagNeumannBoundaryCondition<Type>(bc)
{}

template<class Type>
stagSlipBoundaryCondition<Type>::stagSlipBoundaryCondition
(
    const meshField<Type,staggered>& field,
    const stagSlipBoundaryCondition<Type>& bc
)
:
    stagNeumannBoundaryCondition<Type>(field, bc)
{}

template<class Type>
void stagSlipBoundaryCondition<Type>::eliminateGhosts
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
void stagSlipBoundaryCondition<Type>::eliminateGhosts
(
    linearSystem<symmStencil,Type,staggered>& sys,
    const label l,
    const label d
)
{
    NotImplemented;
}

template<class Type>
void stagSlipBoundaryCondition<Type>::evaluate
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

