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
    const meshField<Type,MeshType>& mshField,
    const boundary& b
)
:
    boundaryCondition<Type,MeshType>(mshField, b),
    neighborProcNum_
    (
        readLabel(b.dict().lookup("neighborProcNum"))
    ),
    tag_
    (
        readLabel(b.dict().lookup("tag"))
    ),
    sendBuffers_(),
    recvBuffers_(),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{
    // Set send/recv buffers for all mesh levels (even though the mesh field may
    // be shallow at this point)

    const labelTensor T(this->T());
    const faceLabel extension(this->extension());

    forAll(this->fvMsh_, l)
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector NSend =
                this->N(l,d)
              + extension.lower()
              + extension.upper();

            const labelVector NRecv = cmptMag(T.T() & NSend);

            sendBuffers_.append(new block<Type>(NSend));
            recvBuffers_.append(new block<Type>(NRecv));
        }
    }
}

template<class Type, class MeshType>
parallelBoundaryCondition<Type,MeshType>::parallelBoundaryCondition
(
    const parallelBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(bc.mshField(), bc.mshBoundary()),
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
    const meshField<Type,MeshType>& field,
    const parallelBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(field, bc.mshBoundary()),
    neighborProcNum_(bc.neighborProcNum_),
    tag_(bc.tag_),
    sendBuffers_(bc.sendBuffers_),
    recvBuffers_(bc.recvBuffers_),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{}

template<class Type, class MeshType>
void parallelBoundaryCondition<Type,MeshType>::initEvaluate
(
    const label l
)
{
    const meshLevel<Type,MeshType>& field = this->mshField()[l];

    const labelVector bo(this->offset());
    const faceLabel extension(this->extension());

    forAll(field, d)
    {
        const meshDirection<Type,MeshType>& fd = field[d];

        const labelVector S(this->S(l,d) - extension.lower());
        const labelVector E(this->E(l,d) + extension.upper());

        block<Type>& sendBuffer =
            sendBuffers_[l*field.size()+d];

        block<Type>& recvBuffer =
            recvBuffers_[l*field.size()+d];

        labelVector ijk;

        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            sendBuffer(ijk-S) = fd(ijk);
        }

        outstandingRecvRequest_ = UPstream::nRequests();

        UIPstream::read
        (
            Pstream::commsTypes::nonBlocking,
            neighborProcNum_,
            reinterpret_cast<char*>(recvBuffer.begin()),
            recvBuffer.byteSize(),
            tag_,
            UPstream::worldComm
        );

        outstandingSendRequest_ = UPstream::nRequests();

        UOPstream::write
        (
            Pstream::commsTypes::nonBlocking,
            neighborProcNum_,
            reinterpret_cast<char*>(sendBuffer.begin()),
            sendBuffer.byteSize(),
            tag_,
            UPstream::worldComm
        );
    }

}

template<class Type, class MeshType>
void parallelBoundaryCondition<Type,MeshType>::evaluate(const label l)
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

    meshLevel<Type,MeshType>& field = this->mshField()[l];

    const labelTensor T(this->T());
    const labelVector bo(this->offset());
    const faceLabel extension(this->extension());

    forAll(field, d)
    {
        meshDirection<Type,MeshType>& fd = field[d];

        const labelVector S(this->S(l,d) - extension.lower());
        const labelVector E(this->E(l,d) + extension.upper());

        block<Type>& recvBuffer = recvBuffers_[l*field.size()+d];

        // Transform the receive buffer back to the orientation of the boundary

        block<Type> data(recvBuffer);
        data.transform(T);

        labelVector ijk;

        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            fd(ijk+bo) = data(ijk-S);
        }
    }
}

}

}

}

