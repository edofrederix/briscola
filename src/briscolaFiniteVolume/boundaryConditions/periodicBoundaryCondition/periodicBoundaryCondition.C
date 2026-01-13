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
    const boundary& b
)
:
    parallelBoundaryCondition<Type,MeshType>(mshField, b)
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
    const meshField<Type,MeshType>& field,
    const periodicBoundaryCondition<Type,MeshType>& bc
)
:
    parallelBoundaryCondition<Type,MeshType>(field, bc)
{}

template<class Type, class MeshType>
void periodicBoundaryCondition<Type,MeshType>::prepare(const label l)
{
    if (this->neighborProcNum_ != Pstream::myProcNo())
        parallelBoundaryCondition<Type,MeshType>::prepare(l);
}

template<class Type, class MeshType>
void periodicBoundaryCondition<Type,MeshType>::evaluate(const label l)
{
    if (this->neighborProcNum_ != Pstream::myProcNo())
    {
        parallelBoundaryCondition<Type,MeshType>::evaluate(l);
    }
    else
    {
        // On the same processor the patch neighbor is assumed to be on the
        // opposing face/edge/vertex. Copy data directly to neighbor. Neighbor
        // will do the same to us.

        meshLevel<Type,MeshType>& field = this->mshField()[l];

        const labelVector bo(this->offset());

        forAll(field, d)
        {
            // Source start and end point

            const labelVector Ss(this->S(l,d));
            const labelVector Es(this->E(l,d));

            // Target start point

            const labelVector St(this->fvMsh_.template S<MeshType>(l,d,-bo));

            labelVector ijk;

            for (ijk.x() = Ss.x(); ijk.x() < Es.x(); ijk.x()++)
            for (ijk.y() = Ss.y(); ijk.y() < Es.y(); ijk.y()++)
            for (ijk.z() = Ss.z(); ijk.z() < Es.z(); ijk.z()++)
            {
                field(d,ijk-Ss+St-bo) = field(d,ijk);
            }
        }
    }
}

}

}

}

