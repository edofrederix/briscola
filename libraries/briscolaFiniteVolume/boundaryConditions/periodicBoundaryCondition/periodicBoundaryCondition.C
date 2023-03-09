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
        // opposing face/edge/vertex. Copy buffers directly.

        const meshField<Type,MeshType>& field = this->mshField();

        forAll(field.boundaryConditions(), bci)
        {
            const boundaryCondition<Type,MeshType>& bc =
                field.boundaryConditions()[bci];

            if (bc.boundaryOffset() == -this->boundaryOffset())
            {
                const labelVector bo(bc.boundaryOffset());
                const meshLevel<Type,MeshType>& fl = field[l];

                forAll(fl, d)
                {
                    const meshDirection<Type,MeshType>& fld = fl[d];

                    labelVector S(fld.boundaryStart(bo));
                    labelVector E(fld.boundaryEnd(bo));

                    block<Type>& recvBuffer = this->recvBuffers_[l*fl.size()+d];

                    labelVector ijk;

                    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                    {
                        recvBuffer(ijk-S) = fld(ijk);
                    }
                }
            }
        }
    }
}

makeBoundaryConditionTypes(periodic)


}

}

}

