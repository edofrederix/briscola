#include "boundaryExchange.H"
#include "periodicBoundary.H"
#include "pMPI.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void boundaryExchange<Type,MeshType>::initBuffers(const label l)
{
    int size = 0;

    labelList neighbors;

    forAll(this->boundaryPtrs_, i)
    {
        const parallelBoundaryCondition<Type,MeshType>& bc =
            *this->boundaryPtrs_[i];

        const faceLabel extension(bc.extension());

        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector S(bc.S(l,d) - extension.lower());
            const labelVector E(bc.E(l,d) + extension.upper());

            size += cmptProduct(E-S);
        }

        neighbors.append(bc.mshBoundary().neighborProcNum());
    }

    sendBuffers_.set(l, new List<Type>(size));
    recvBuffers_.set(l, new List<Type>(size));
}

template<class Type, class MeshType>
void boundaryExchange<Type,MeshType>::correct
(
    const label l,
    const label baseType
)
{
    if (!sendBuffers_.set(l))
        this->initBuffers(l);

    List<Type>& sendBuffer = sendBuffers_[l];
    List<Type>& recvBuffer = recvBuffers_[l];

    meshLevel<Type,MeshType>& field = field_[l];

    labelList counts(this->boundaryPtrs_.size(), 0);
    labelList displacements(this->boundaryPtrs_.size());

    // Copy data to send buffer

    int c = 0;
    forAll(this->boundaryPtrs_, i)
    {
        const parallelBoundaryCondition<Type,MeshType>& bc =
            *this->boundaryPtrs_[i];

        displacements[i] = c;

        if (baseType == -1 || bc.baseType() == baseType)
        {
            const faceLabel extension(bc.extension());

            for (int d = 0; d < MeshType::numberOfDirections; d++)
            {
                const labelVector S(bc.S(l,d) - extension.lower());
                const labelVector E(bc.E(l,d) + extension.upper());

                labelVector ijk;
                for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                {
                    sendBuffer[c++] = field(ijk);
                }

                counts[i] += cmptProduct(E-S);
            }
        }
    }

    // Exchange data

    MPI_Neighbor_alltoallv
    (
        sendBuffer.begin(),
        counts.begin(),
        displacements.begin(),
        mpiDataTypePtr_->type,
        recvBuffer.begin(),
        counts.begin(),
        displacements.begin(),
        mpiDataTypePtr_->type,
        fvMsh_.msh().comm()
    );

    // Unpack data

    c = 0;
    forAll(this->boundaryPtrs_, i)
    {
        const parallelBoundaryCondition<Type,MeshType>& bc =
            *this->boundaryPtrs_[i];

        if (baseType == -1 || bc.baseType() == baseType)
        {
            const faceLabel extension(bc.extension());
            const labelVector bo(bc.offset());

            for (int d = 0; d < MeshType::numberOfDirections; d++)
            {
                const labelVector S(bc.S(l,d) - extension.lower());
                const labelVector E(bc.E(l,d) + extension.upper());

                labelVector ijk;
                for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                {
                    field(ijk+bo) = recvBuffer[c++];
                }
            }
        }
    }
}

template<class Type, class MeshType>
boundaryExchange<Type,MeshType>::boundaryExchange
(
    meshField<Type,MeshType>& field
)
:
    field_(field),
    fvMsh_(field.fvMsh_),
    sendBuffers_(fvMsh_.size()),
    recvBuffers_(fvMsh_.size()),
    mpiDataTypePtr_(new pMPI<Type>())
{
    typedef parallelBoundaryCondition<Type,MeshType> pbcType;

    forAll(field_.boundaryConditions(), i)
        if (field_.boundaryConditions()[i].template castable<pbcType>())
            boundaryPtrs_.append
            (
                &field_.boundaryConditions()[i].template cast<pbcType>()
            );
}

template<class Type, class MeshType>
boundaryExchange<Type,MeshType>::boundaryExchange
(
    boundaryExchange<Type,MeshType>& s
)
:
    field_(s.field_),
    fvMsh_(s.fvMsh_),
    boundaryPtrs_(s.boundaryPtrs_),
    sendBuffers_(s.sendBuffers_),
    recvBuffers_(s.recvBuffers_),
    mpiDataTypePtr_(new pMPI<Type>())
{}

template<class Type, class MeshType>
boundaryExchange<Type,MeshType>::~boundaryExchange()
{}

// Instantiate

makeboundaryExchanges(colocated)
makeboundaryExchanges(staggered)

}

}

}
