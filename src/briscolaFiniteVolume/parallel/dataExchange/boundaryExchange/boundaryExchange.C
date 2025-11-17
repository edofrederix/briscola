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
    int size = 0;

    // Allocate for all levels and directions. Buffers may only be used
    // partially, depending on which levels are corrected. The overhead of
    // having multiple levels is relatively small.

    forAll(this->boundaryPtrs_, i)
    {
        const parallelBoundaryCondition<Type,MeshType>& bc =
            *this->boundaryPtrs_[i];

        for (int l = 0; l < fvMsh_.size(); l++)
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const labelVector S(bc.S(l,d));
            const labelVector E(bc.E(l,d));

            size += cmptProduct(E-S);
        }
    }

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
    if (!sendBuffers_.size())
        this->initBuffers();

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
            forAll(field_, L)
            {
                // Only for the selected level, unless all levels are selected
                // (negative l value).

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

                        counts[i] += cmptProduct(E-S);
                    }
                }
            }
        }
    }

    // Exchange data

    if (Pstream::parRun())
    {
        MPI_Neighbor_alltoallv
        (
            sendBuffers_.begin(),
            counts.begin(),
            displacements.begin(),
            mpiDataTypePtr_->type,
            recvBuffers_.begin(),
            counts.begin(),
            displacements.begin(),
            mpiDataTypePtr_->type,
            PstreamGlobals::MPI_COMM_FOAM
        );
    }
    else
    {
        // When running in serial, the only exchange can occur across periodic
        // boundaries. Copy manually.

        if (baseType == -1 || baseType == PERIODICBC)
        forAll(this->boundaryPtrs_, i)
        {
            const parallelBoundaryCondition<Type,MeshType>& bc1 =
                *this->boundaryPtrs_[i];

            const labelVector bo1(bc1.offset());

            forAll(this->boundaryPtrs_, j)
            {
                const parallelBoundaryCondition<Type,MeshType>& bc2 =
                    *this->boundaryPtrs_[j];

                const labelVector bo2(bc2.offset());

                // Find opposite boundary and copy

                if (bo1 == -bo2)
                {
                    for (int k = 0; k < counts[i]; k++)
                        recvBuffers_[displacements[j]+k] =
                            sendBuffers_[displacements[i]+k];

                    break;
                }
            }
        }
    }

    // Unpack data

    c = 0;
    forAll(this->boundaryPtrs_, i)
    {
        const parallelBoundaryCondition<Type,MeshType>& bc =
            *this->boundaryPtrs_[i];

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

                            forAllBlockLinear(B, j)
                                B(j) = recvBuffers_[c++];

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
    recvBuffers_(0),
    mpiDataTypePtr_(new pMPI<Type>())
{
    typedef parallelBoundaryCondition<Type,MeshType> pbcType;

    if (!field_.boundaryConditions().size())
        field_.addBoundaryConditions();

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
