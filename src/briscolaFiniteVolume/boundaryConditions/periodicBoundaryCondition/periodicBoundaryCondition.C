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
    const meshLevel<Type,MeshType>& level,
    const boundary& b
)
:
    parallelBoundaryCondition<Type,MeshType>(level, b)
{}

template<class Type, class MeshType>
periodicBoundaryCondition<Type,MeshType>::periodicBoundaryCondition
(
    const periodicBoundaryCondition<Type,MeshType>& bc
)
:
    parallelBoundaryCondition<Type,MeshType>(bc)
{}

template<class Type, class MeshType>
periodicBoundaryCondition<Type,MeshType>::periodicBoundaryCondition
(
    const periodicBoundaryCondition<Type,MeshType>& bc,
    const meshLevel<Type,MeshType>& level
)
:
    parallelBoundaryCondition<Type,MeshType>(bc, level)
{}

template<class Type, class MeshType>
void periodicBoundaryCondition<Type,MeshType>::prepare()
{
    if (this->neighborProcNum_ != Pstream::myProcNo())
        parallelBoundaryCondition<Type,MeshType>::prepare();
}

template<class Type, class MeshType>
void periodicBoundaryCondition<Type,MeshType>::evaluate()
{
    if (this->neighborProcNum_ != Pstream::myProcNo())
    {
        parallelBoundaryCondition<Type,MeshType>::evaluate();
    }
    else
    {
        // On the same processor the patch neighbor is assumed to be on the
        // opposing face/edge/vertex. Copy data directly to neighbor. Neighbor
        // will do the same to us.

        const labelVector bo(this->offset());
        const label l(this->l_);

        meshLevel<Type,MeshType>& level = this->level_;

        forAll(level, d)
        {
            // Source start and end point

            const labelVector Ss(this->S(d));
            const labelVector Es(this->E(d));

            // Target start point

            const labelVector St(this->fvMsh_.template S<MeshType>(l,d,-bo));

            labelVector ijk;

            for (ijk.x() = Ss.x(); ijk.x() < Es.x(); ijk.x()++)
            for (ijk.y() = Ss.y(); ijk.y() < Es.y(); ijk.y()++)
            for (ijk.z() = Ss.z(); ijk.z() < Es.z(); ijk.z()++)
            {
                level(d,ijk-Ss+St-bo) = level(d,ijk);
            }
        }
    }
}

}

}

}

