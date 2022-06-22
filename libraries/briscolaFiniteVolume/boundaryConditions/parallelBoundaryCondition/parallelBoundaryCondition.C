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
    const partPatch& patch
)
:
    boundaryCondition<Type,MeshType>(mshField, patch),
    neighborProcNum_
    (
        readLabel(patch.dict().lookup("neighborProcNum"))
    ),
    tag_
    (
        readLabel(patch.dict().lookup("tag"))
    ),
    sendBuffers_(),
    recvBuffers_(),
    order_(MeshType::numberOfDirections, -1),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)

{
    const labelVector bo(this->boundaryOffset());

    // Each direction will be received in the order of directions of the
    // neighbor. We should process data in that order. A direction is uniquely
    // identified by its padding, so we find the mapping between local direction
    // (d1) and neighbor direction (d2) by transforming paddings.

    const labelTensor T(this->T());

    forAll(mshField[0], d2)
    {
        const labelVector padding
        (
            Foam::cmptMag(T & MeshType::padding[d2])
        );

        forAll(mshField[0], d1)
        {
            if (padding == MeshType::padding[d1])
            {
                order_[d2] = d1;
                break;
            }
        }

        if (order_[d2] == -1)
        {
            FatalErrorInFunction
                << "Could not determine mesh direction ordering at neighbor" << endl
                << exit(FatalError);
        }
    }

    forAll(mshField, l)
    {
        forAll(mshField[l], d)
        {
            const label d1 = order_[d];

            const labelVector NSend =
                mshField[l][d].extendedBoundaryN(bo);

            const labelVector NRecv =
                cmptMag(T.T() & mshField[l][d1].extendedBoundaryN(bo));

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
    boundaryCondition<Type,MeshType>(bc.mshField(), bc.patch()),
    neighborProcNum_(bc.neighborProcNum_),
    tag_(bc.tag_),
    sendBuffers_(bc.sendBuffers_),
    recvBuffers_(bc.recvBuffers_),
    order_(bc.order_),
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
    boundaryCondition<Type,MeshType>(field, bc.patch()),
    neighborProcNum_(bc.neighborProcNum_),
    tag_(bc.tag_),
    sendBuffers_(bc.sendBuffers_),
    recvBuffers_(bc.recvBuffers_),
    order_(bc.order_),
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

    const labelVector bo(this->boundaryOffset());

    forAll(field, d)
    {
        const meshDirection<Type,MeshType>& fd = field[d];

        const labelVector S(fd.extendedBoundaryStart(bo));
        const labelVector E(fd.extendedBoundaryEnd(bo));
        const labelVector C(fd.copyOffset(bo));

        block<Type>& sendBuffer =
            sendBuffers_[l*field.size()+d];

        block<Type>& recvBuffer =
            recvBuffers_[l*field.size()+d];

        labelVector ijk;

        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            sendBuffer(ijk-S) = fd(ijk+C);
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
void parallelBoundaryCondition<Type,MeshType>::evaluate
(
    const label l,
    const bool
)
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
    const labelVector bo(this->boundaryOffset());

    forAll(field, d2)
    {
        // Direction d1 corresponds to the neighbor's direction d2

        const label d1 = order_[d2];

        meshDirection<Type,MeshType>& fd = field[d1];

        const labelVector S(fd.extendedBoundaryStart(bo));
        const labelVector E(fd.extendedBoundaryEnd(bo));

        block<Type>& recvBuffer =
            recvBuffers_[l*field.size()+d2];

        // Transform the receive buffer back to the orientation of the patch

        block<Type> data(recvBuffer);
        data.transform(T);

        // For parallel on-boundary BCs, the linear system is not solved
        // and the solution set to the given value. For off-boundary BCs,
        // set the ghost cell appropriately

        labelVector ijk;

        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            fd(ijk+bo) = data(ijk-S);
        }
    }
}

makeBoundaryConditionType(parallel,label,colocated)
makeBoundaryConditionType(parallel,label,staggered)

makeBoundaryConditionType(parallel,scalar,colocated)
makeBoundaryConditionType(parallel,scalar,staggered)

makeBoundaryConditionType(parallel,hexScalar,colocated)
makeBoundaryConditionType(parallel,hexScalar,staggered)

makeBoundaryConditionType(parallel,vector,colocated)
makeBoundaryConditionType(parallel,vector,staggered)

makeBoundaryConditionType(parallel,hexVector,colocated)
makeBoundaryConditionType(parallel,hexVector,staggered)

makeBoundaryConditionType(parallel,tensor,colocated)
makeBoundaryConditionType(parallel,tensor,staggered)

makeBoundaryConditionType(parallel,sphericalTensor,colocated)
makeBoundaryConditionType(parallel,sphericalTensor,staggered)

makeBoundaryConditionType(parallel,symmTensor,colocated)
makeBoundaryConditionType(parallel,symmTensor,staggered)

makeBoundaryConditionType(parallel,diagTensor,colocated)
makeBoundaryConditionType(parallel,diagTensor,staggered)

makeBoundaryConditionType(parallel,stencil,colocated)
makeBoundaryConditionType(parallel,stencil,staggered)

makeBoundaryConditionType(parallel,symmStencil,colocated)
makeBoundaryConditionType(parallel,symmStencil,staggered)

makeBoundaryConditionType(parallel,diagStencil,colocated)
makeBoundaryConditionType(parallel,diagStencil,staggered)

}

}

}

