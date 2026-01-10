#include "parallelBoundaryCondition.H"

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
parallelBoundaryCondition<Type,MeshType>::parallelBoundaryCondition
(
    const meshLevel<Type,MeshType>& level,
    const boundary& b
)
:
    boundaryCondition<Type,MeshType>(level, b),
    neighborProcNum_
    (
        b.cast<parallelBoundary>().neighborProcNum()
    ),
    tag_
    (
        b.cast<parallelBoundary>().tag()
    ),
    sendBuffers_(MeshType::numberOfDirections),
    recvBuffers_(MeshType::numberOfDirections),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{
    // Set send/recv buffers

    const labelTensor T(this->T());

    for (int d = 0; d < MeshType::numberOfDirections; d++)
    {
        const labelVector NSend = this->N(d);
        const labelVector NRecv = cmptMag(T.T() & NSend);

        sendBuffers_.set(d, new block<Type>(NSend));
        recvBuffers_.set(d, new block<Type>(NRecv));
    }
}

template<class Type, class MeshType>
parallelBoundaryCondition<Type,MeshType>::parallelBoundaryCondition
(
    const parallelBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(bc),
    neighborProcNum_(bc.neighborProcNum_),
    tag_(bc.tag_),
    sendBuffers_(bc.sendBuffers_),
    recvBuffers_(bc.recvBuffers_),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{}

template<class Type, class MeshType>
parallelBoundaryCondition<Type,MeshType>::parallelBoundaryCondition
(
    const meshLevel<Type,MeshType>& level,
    const parallelBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(level, bc),
    neighborProcNum_(bc.neighborProcNum_),
    tag_(bc.tag_),
    sendBuffers_(bc.sendBuffers_),
    recvBuffers_(bc.recvBuffers_),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{}

template<class Type, class MeshType>
void parallelBoundaryCondition<Type,MeshType>::prepare()
{
    const labelVector bo(this->offset());

    const meshLevel<Type,MeshType>& level = this->level_;

    forAll(level, d)
    {
        const labelVector S(this->S(d));
        const labelVector E(this->E(d));

        block<Type>& sendBuffer = sendBuffers_[d];
        block<Type>& recvBuffer = recvBuffers_[d];

        labelVector ijk;

        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            sendBuffer(ijk-S) = level(d,ijk);
        }

        outstandingRecvRequest_ = UPstream::nRequests();

        UIPstream::read
        (
            Pstream::commsTypes::nonBlocking,
            neighborProcNum_,
            reinterpret_cast<char*>(recvBuffer.begin()),
            recvBuffer.byteSize(),
            tag_,
            level.lvl().comms()
        );

        outstandingSendRequest_ = UPstream::nRequests();

        UOPstream::write
        (
            Pstream::commsTypes::nonBlocking,
            neighborProcNum_,
            reinterpret_cast<char*>(sendBuffer.begin()),
            sendBuffer.byteSize(),
            tag_,
            level.lvl().comms()
        );
    }
}

template<class Type, class MeshType>
void parallelBoundaryCondition<Type,MeshType>::evaluate()
{
    if
    (
        outstandingRecvRequest_ >= 0
     && outstandingRecvRequest_ < Pstream::nRequests()
    )
    {
        UPstream::waitRequest(outstandingRecvRequest_);
    }

    outstandingSendRequest_ = -1;
    outstandingRecvRequest_ = -1;

    const labelTensor T(this->T());
    const labelVector bo(this->offset());

    meshLevel<Type,MeshType>& level = this->level_;

    forAll(level, d)
    {
        const labelVector S(this->S(d));
        const labelVector E(this->E(d));

        block<Type>& recvBuffer = recvBuffers_[d];

        // Transform the receive buffer back to the orientation of the boundary

        block<Type> data(recvBuffer);
        data.transform(T);

        labelVector ijk;

        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            level(d,ijk+bo) = data(ijk-S);
        }
    }
}

}

}

}

