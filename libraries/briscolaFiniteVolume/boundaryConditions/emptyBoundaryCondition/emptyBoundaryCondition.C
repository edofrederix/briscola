#include "emptyBoundaryCondition.H"

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
emptyBoundaryCondition<Type,MeshType>::emptyBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b
)
:
    boundaryCondition<Type,MeshType>(mshField, b)
{
    // Empty boundaries can only sit normal to unit mesh size. Check.

    const part& p = this->fvMsh_[0];
    const labelVector bo(this->offset());

    for (label d = 0; d < 3; d++)
        if (Foam::mag(bo[d]) != 0 && p.N()[d] != 1)
            FatalErrorInFunction
                << "The empty boundary condition can only be applied "
                << "normal to mesh dimensions with only 1 cell thickness."
                << endl << exit(FatalError);
}

template<class Type, class MeshType>
emptyBoundaryCondition<Type,MeshType>::emptyBoundaryCondition
(
    const emptyBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(bc.mshField(), bc.mshBoundary())
{}

template<class Type, class MeshType>
emptyBoundaryCondition<Type,MeshType>::emptyBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const emptyBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(field, bc.mshBoundary())
{}

template<class Type, class MeshType>
void emptyBoundaryCondition<Type,MeshType>::prepare(const label)
{}

template<class Type, class MeshType>
void emptyBoundaryCondition<Type,MeshType>::evaluate(const label l)
{
    meshLevel<Type,MeshType>& field = this->mshField()[l];

    const labelVector bo(this->offset());
    const faceLabel extension(this->extension());

    forAll(field, d)
    {
        meshDirection<Type,MeshType>& fd = field[d];

        const labelVector S(this->S(l,d) - extension.lower());
        const labelVector E(this->E(l,d) + extension.upper());

        // For shifted boundaries, the boundary values are set to zero. For
        // non-shifted boundaries, apply homogeneous Neumann

        labelVector ijk;

        if (MeshType::shifted(d,bo))
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                fd(ijk+bo) = Zero;
            }
        }
        else
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                fd(ijk+bo) = fd(ijk);
            }
        }
    }
}

}

}

}

