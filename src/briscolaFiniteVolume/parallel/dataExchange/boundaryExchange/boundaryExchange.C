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

    // Allocate for all directions

    int size = 0;

    forAll(this->boundaryPtrs_, i)
    if (this->boundaryPtrs_[i])
    {
        const pbcType& bc = *this->boundaryPtrs_[i];

        if (bc.neighborProcNum() != Pstream::myProcNo())
            for (int d = 0; d < MeshType::numberOfDirections; d++)
                size += cmptProduct(bc.E(d) - bc.S(d));
    }

    // Exactly the same amount of data is exchanged so send and receive buffers
    // are the same length

    sendBuffers_.resize(size);
    recvBuffers_.resize(size);
}

template<class Type, class MeshType>
void boundaryExchange<Type,MeshType>::correct(const label baseType)
{
    typedef parallelBoundaryCondition<Type,MeshType> pbcType;

    const levelComms& comms = level_.lvl().comms();

    if (!sendBuffers_.size())
        this->initBuffers();

    // Set send data

    labelList sendCount(comms.sendNeighbors().size(), 0);
    labelList sendDisp(comms.sendNeighbors().size(), 0);

    label c = 0;
    forAll(comms.sendMap(), i)
    {
        forAll(comms.sendMap()[i], j)
        {
            const label k = comms.sendMap()[i][j];
            const pbcType& bc = *this->boundaryPtrs_[k];

            if (baseType == -1 || bc.baseType() == baseType)
            {
                for (int d = 0; d < MeshType::numberOfDirections; d++)
                {
                    const labelVector S(bc.S(d));
                    const labelVector E(bc.E(d));

                    labelVector ijk;
                    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                    {
                        sendBuffers_[c++] = level_(d,ijk);
                    }

                    sendCount[i] += cmptProduct(E-S);
                }
            }
        }

        if (i > 0)
            sendDisp[i] = sendDisp[i-1] + sendCount[i-1];
    }

    // Set receive data

    labelList recvCount(comms.recvNeighbors().size(), 0);
    labelList recvDisp(comms.recvNeighbors().size(), 0);

    forAll(comms.recvMap(), i)
    {
        forAll(comms.recvMap()[i], j)
        {
            const label k = comms.recvMap()[i][j];
            const pbcType& bc = *this->boundaryPtrs_[k];

            if (baseType == -1 || bc.baseType() == baseType)
                for (int d = 0; d < MeshType::numberOfDirections; d++)
                    recvCount[i] +=
                        cmptProduct(bc.E(d) - bc.S(d));
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
            comms.communicator()
        );
    }

    // Unpack data

    c = 0;
    forAll(comms.recvMap(), i)
    forAll(comms.recvMap()[i], j)
    {
        const label k = comms.recvMap()[i][j];
        const pbcType& bc = *this->boundaryPtrs_[k];

        if (baseType == -1 || bc.baseType() == baseType)
        {
            const labelVector bo(bc.offset());
            const labelTensor T(bc.T());

            if (T == eye)
            {
                for (int d = 0; d < MeshType::numberOfDirections; d++)
                {
                    const labelVector S(bc.S(d));
                    const labelVector E(bc.E(d));

                    labelVector ijk;
                    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                    {
                        level_(d,ijk+bo) = recvBuffers_[c++];
                    }
                }
            }
            else
            {
                // Store in intermediate block to transform the data

                for (int d = 0; d < MeshType::numberOfDirections; d++)
                {
                    const labelVector S(bc.S(d));
                    const labelVector E(bc.E(d));
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
                        level_(d,ijk+bo) = B(ijk-S);
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
                        this->boundaryPtrs_[i]->evaluate();
}

template<class Type, class MeshType>
boundaryExchange<Type,MeshType>::boundaryExchange
(
    meshLevel<Type,MeshType>& level
)
:
    level_(level),
    fvMsh_(level.fvMsh_),
    sendBuffers_(0),
    recvBuffers_(0)
{
    if (Pstream::parRun())
        mpiDataTypePtr_.reset(new pMPI<Type>());

    typedef parallelBoundaryCondition<Type,MeshType> pbcType;

    if (!level_.boundaryConditions().size())
        level_.addBoundaryConditions();

    boundaryPtrs_.resize(level_.boundaryConditions().size());

    forAll(level_.boundaryConditions(), i)
    {
        if (level_.boundaryConditions()[i].template castable<pbcType>())
        {
            boundaryPtrs_[i] =
                &level_.boundaryConditions()[i].template cast<pbcType>();
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
    level_(s.level_),
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
