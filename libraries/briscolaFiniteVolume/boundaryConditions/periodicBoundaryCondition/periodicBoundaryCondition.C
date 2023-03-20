#include "periodicBoundaryCondition.H"

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
periodicBoundaryCondition<Type,MeshType>::periodicBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const partPatch& patch
)
:
    parallelBoundaryCondition<Type,MeshType>(mshField, patch)
{}

template<class Type, class MeshType>
periodicBoundaryCondition<Type,MeshType>::periodicBoundaryCondition
(
    const periodicBoundaryCondition<Type,MeshType>& bc
)
:
    parallelBoundaryCondition<Type,MeshType>(bc.mshField(), bc.patch())
{}

template<class Type, class MeshType>
periodicBoundaryCondition<Type,MeshType>::periodicBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const periodicBoundaryCondition<Type,MeshType>& bc
)
:
    parallelBoundaryCondition<Type,MeshType>(field, bc.patch())
{}

template<class Type, class MeshType>
void periodicBoundaryCondition<Type,MeshType>::initEvaluate
(
    const label l
)
{
    if (this->neighborProcNum_ != Pstream::myProcNo())
    {
        parallelBoundaryCondition<Type,MeshType>::initEvaluate(l);
    }
    else
    {
        // On the same processor the patch neighbor is assumed to be on the
        // opposing face/edge/vertex. Copy data directly to neighbor. Neighbor
        // will do the same to us.

        meshLevel<Type,MeshType>& field = this->mshField()[l];

        const labelVector bo(this->boundaryOffset());

        forAll(field, d)
        {
            meshDirection<Type,MeshType>& fd = field[d];

            // Source and target start point

            const labelVector Ss
            (
                fd.activeBoundaryStart(bo)
              - this->extension_.lower()
            );

            const labelVector St
            (
                fd.activeBoundaryStart(-bo)
              - this->extension_.lower()
            );

            const labelVector Es
            (
                fd.activeBoundaryEnd(bo)
              + this->extension_.upper()
            );

            labelVector ijk;

            for (ijk.x() = Ss.x(); ijk.x() < Es.x(); ijk.x()++)
            for (ijk.y() = Ss.y(); ijk.y() < Es.y(); ijk.y()++)
            for (ijk.z() = Ss.z(); ijk.z() < Es.z(); ijk.z()++)
            {
                fd(ijk-Ss+St-bo) = fd(ijk);
            }
        }
    }
}

template<class Type, class MeshType>
void periodicBoundaryCondition<Type,MeshType>::evaluate
(
    const label l,
    const bool
)
{
    if (this->neighborProcNum_ != Pstream::myProcNo())
    {
        parallelBoundaryCondition<Type,MeshType>::evaluate(l);
    }
}

}

}

}

