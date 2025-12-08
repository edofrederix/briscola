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
void boundaryExchange<Type,MeshType>::initBuffers()
{
    typedef parallelBoundaryCondition<Type,MeshType> pbcType;

    // Allocate for all levels and directions. Buffers may only be used
    // partially, depending on which levels are corrected. The overhead of
    // having multiple levels is relatively small.

    int size = 0;

    forAll(this->boundaryPtrs_, i)
    if (this->boundaryPtrs_[i])
    {
        const pbcType& bc = *this->boundaryPtrs_[i];

        if (bc.neighborProcNum() != Pstream::myProcNo())
            for (int l = 0; l < fvMsh_.size(); l++)
                for (int d = 0; d < MeshType::numberOfDirections; d++)
                    size += cmptProduct(bc.E(l,d) - bc.S(l,d));
    }

    // Exactly the same amount of data is exchanged so send and receive buffers
    // are the same length

    sendBuffers_.resize(size);
    recvBuffers_.resize(size);
}

template<class Type, class MeshType>
void boundaryExchange<Type,MeshType>::correct
(
    const label l,
    const label baseType
)
{
    typedef parallelBoundaryCondition<Type,MeshType> pbcType;

    if (!sendBuffers_.size())
        this->initBuffers();

    // Set send data

    labelList sendCount(fvMesh::MPI_neighbors_send.size(), 0);
    labelList sendDisp(fvMesh::MPI_neighbors_send.size(), 0);

    label c = 0;
    forAll(fvMesh::MPI_neighbors_send_map, i)
    {
        forAll(fvMesh::MPI_neighbors_send_map[i], j)
        {
            const label k = fvMesh::MPI_neighbors_send_map[i][j];
            const pbcType& bc = *this->boundaryPtrs_[k];

            if (baseType == -1 || bc.baseType() == baseType)
            {
                forAll(field_, L)
                {
                    // Only for the selected level, unless all levels are
                    // selected (negative l value).

                    if (L == l || l < 0)
                    {
                        for (int d = 0; d < MeshType::numberOfDirections; d++)
                        {
                            const labelVector S(bc.S(L,d));
                            const labelVector E(bc.E(L,d));

                            labelVector ijk;
                            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                            {
                                sendBuffers_[c++] = field_[L](d,ijk);
                            }

                            sendCount[i] += cmptProduct(E-S);
                        }
                    }
                }
            }
        }

        if (i > 0)
            sendDisp[i] = sendDisp[i-1] + sendCount[i-1];
    }

    // Set receive data

    labelList recvCount(fvMesh::MPI_neighbors_recv.size(), 0);
    labelList recvDisp(fvMesh::MPI_neighbors_recv.size(), 0);

    forAll(fvMesh::MPI_neighbors_recv_map, i)
    {
        forAll(fvMesh::MPI_neighbors_recv_map[i], j)
        {
            const label k = fvMesh::MPI_neighbors_recv_map[i][j];
            const pbcType& bc = *this->boundaryPtrs_[k];

            if (baseType == -1 || bc.baseType() == baseType)
                forAll(field_, L)
                    if (L == l || l < 0)
                        for (int d = 0; d < MeshType::numberOfDirections; d++)
                            recvCount[i] +=
                                cmptProduct(bc.E(L,d) - bc.S(L,d));
        }

        if (i > 0)
            recvDisp[i] = recvDisp[i-1] + recvCount[i-1];
    }

    // Exchange data

    if (Pstream::parRun())
    {
        MPI_Neighbor_alltoallv
        (
            sendBuffers_.begin(),
            sendCount.begin(),
            sendDisp.begin(),
            mpiDataTypePtr_->type,
            recvBuffers_.begin(),
            recvCount.begin(),
            recvDisp.begin(),
            mpiDataTypePtr_->type,
            PstreamGlobals::MPI_COMM_FOAM
        );
    }

    // Unpack data

    c = 0;
    forAll(fvMesh::MPI_neighbors_recv_map, i)
    forAll(fvMesh::MPI_neighbors_recv_map[i], j)
    {
        const label k = fvMesh::MPI_neighbors_recv_map[i][j];
        const pbcType& bc = *this->boundaryPtrs_[k];

        if (baseType == -1 || bc.baseType() == baseType)
        {
            const labelVector bo(bc.offset());
            const labelTensor T(bc.T());

            forAll(field_, L)
            {
                // Only for the selected level, unless all levels are selected
                // (negative l value).

                if (L == l || l < 0)
                {
                    if (T == eye)
                    {
                        for (int d = 0; d < MeshType::numberOfDirections; d++)
                        {
                            const labelVector S(bc.S(L,d));
                            const labelVector E(bc.E(L,d));

                            labelVector ijk;
                            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                            {
                                field_[L](d,ijk+bo) = recvBuffers_[c++];
                            }
                        }
                    }
                    else
                    {
                        // Store in intermediate block to transform the data

                        for (int d = 0; d < MeshType::numberOfDirections; d++)
                        {
                            const labelVector S(bc.S(L,d));
                            const labelVector E(bc.E(L,d));
                            const labelVector N(cmptMag(T.T() & (E-S)));

                            block<Type> B(N);

                            forAllBlockLinear(B, ii)
                                B(ii) = recvBuffers_[c++];

                            B.transform(T);

                            labelVector ijk;
                            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                            {
                                field_[L](d,ijk+bo) = B(ijk-S);
                            }
                        }
                    }
                }
            }
        }
    }

    // Handle periodic boundaries to self

    const label myProcNo = Pstream::myProcNo();

    forAll(this->boundaryPtrs_, i)
        if (this->boundaryPtrs_[i])
            if (this->boundaryPtrs_[i]->neighborProcNum() == myProcNo)
                if (baseType == -1 || baseType == PERIODICBC)
                    forAll(field_, L)
                        if (L == l || l < 0)
                            this->boundaryPtrs_[i]->evaluate(L);
}

template<class Type, class MeshType>
boundaryExchange<Type,MeshType>::boundaryExchange
(
    meshField<Type,MeshType>& field
)
:
    field_(field),
    fvMsh_(field.fvMsh_),
    sendBuffers_(0),
    recvBuffers_(0)
{
    if (Pstream::parRun())
        mpiDataTypePtr_.reset(new pMPI<Type>());

    typedef parallelBoundaryCondition<Type,MeshType> pbcType;

    if (!field_.boundaryConditions().size())
        field_.addBoundaryConditions();

    boundaryPtrs_.resize(field_.boundaryConditions().size());

    forAll(field_.boundaryConditions(), i)
    {
        if (field_.boundaryConditions()[i].template castable<pbcType>())
        {
            boundaryPtrs_[i] =
                &field_.boundaryConditions()[i].template cast<pbcType>();
        }
        else
        {
            boundaryPtrs_[i] = nullptr;
        }
    }
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
    recvBuffers_(s.recvBuffers_)
{
    if (Pstream::parRun())
        mpiDataTypePtr_.reset(new pMPI<Type>());
}

template<class Type, class MeshType>
boundaryExchange<Type,MeshType>::~boundaryExchange()
{}

// Instantiate

makeboundaryExchanges(colocated)
makeboundaryExchanges(staggered)

}

}

}
