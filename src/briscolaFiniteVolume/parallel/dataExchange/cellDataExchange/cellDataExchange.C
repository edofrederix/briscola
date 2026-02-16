#include "cellDataExchange.H"
#include "patchBoundary.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTemplateTypeNameAndDebug(cellDataExchange<colocated>, 0);
defineTemplateTypeNameAndDebug(cellDataExchange<staggered>, 0);

template<class MeshType>
cellDataExchange<MeshType>::cellDataExchange
(
    const List<labelVector>& indices,
    const fvMesh& fvMsh,
    const label l,
    const label d
)
:
    dataExchange<MeshType>(fvMsh,l,d)
{
    init(indices);
}

template<class MeshType>
cellDataExchange<MeshType>::cellDataExchange
(
    const cellDataExchange& e
)
:
    dataExchange<MeshType>(e),
    indices_(e.indices_),
    map_(e.map_),
    neighbors_(e.neighbors_),
    neighborOffsets_(e.neighborOffsets_),
    cellsToSend_(e.cellsToSend_),
    cellsToRecv_(e.cellsToRecv_),
    sendTags_(e.sendTags_)
{}

template<class MeshType>
cellDataExchange<MeshType>::~cellDataExchange()
{}

template<class MeshType>
void cellDataExchange<MeshType>::init(const List<labelVector>& indices)
{
    indices_ = indices;

    const faceLabel I = this->fvMsh_.template I<MeshType>(this->l_,this->d_);

    const labelVector U = I.upper();
    const labelVector L = I.lower();

    const level& lvl = this->fvMsh_[this->l_];

    // Set neighboring processor numbers, offsets and transformations. The first
    // one is us.

    neighbors_.clear();
    neighbors_.setSize(lvl.boundaries().size()+1,-1);
    neighbors_[0] = Pstream::myProcNo();

    neighborOffsets_.clear();
    neighborOffsets_.setSize(neighbors_.size(), zeroXYZ);

    List<labelTensor> neighborTs(neighbors_.size(), eye);

    forAll(lvl.boundaries(), i)
    {
        if (lvl.boundaries()[i].castable<parallelBoundary>())
        {
            const parallelBoundary& b =
                lvl.boundaries()[i].cast<parallelBoundary>();

            neighbors_[i+1] = b.neighborProcNum();
            neighborOffsets_[i+1] = b.offset();
            neighborTs[i+1] = b.T();
        }
    }

    // Collect the required cells and store per neighbor processor

    cellsToRecv_.clear();
    cellsToRecv_.setSize(lvl.boundaries().size()+1);

    map_.clear();
    map_.setSize(indices.size());

    forAll(indices, i)
    {
        labelVector index = indices[i];

        const labelVector closest = labelVector
        (
            Foam::max(Foam::min(index.x(), U.x()), L.x()-1),
            Foam::max(Foam::min(index.y(), U.y()), L.y()-1),
            Foam::max(Foam::min(index.z(), U.z()), L.z()-1)
        );

        const labelVector offset = index - closest;

        if (offset == zeroXYZ)
        {
            // Local point

            cellsToRecv_[0].append(index);
            map_[i] = Tuple2<label,label>(0,cellsToRecv_[0].size()-1);
        }
        else
        {
            const labelVector bo =
                briscola::cmptMin(briscola::cmptMax(offset, -unitXYZ), unitXYZ);

            const label degree = cmptSum(cmptMag(bo));

            if (!this->fvMsh_.structured() && degree > 1)
                FatalErrorInFunction
                    << "Can only obtain cell data across faces "
                    << "on unstructured meshes" << endl
                    << abort(FatalError);

            // Get the boundary to the other processor

            const label bNum = findIndex(neighborOffsets_, bo);

            // Error when not found

            if (bNum == -1)
                FatalErrorInFunction
                    << "Could not find boundary corresponding to offset "
                    << bo << endl << abort(FatalError);

            // Error when we're outside the domain

            if
            (
                lvl.boundaries()[bNum-1].castable<patchBoundary>()
             || lvl.boundaries()[bNum-1].castable<emptyBoundary>()
            )
                FatalErrorInFunction
                    << "Cannot exchange cell data across a domain boundary "
                    << "at offset " << bo << endl << abort(FatalError);

            // Store trimmed index

            for (int j = 0; j < 3; j++)
                if (bo[j] != 0)
                    index[j] -= bo[j] > 0 ? (U[j]-1) :L[j];

            cellsToRecv_[bNum].append(index);

            map_[i] =
                Tuple2<label,label>
                (
                    bNum,
                    cellsToRecv_[bNum].size()-1
                );
        }
    }

    // Exchange processor numbers, sizes, boundary offsets and tags using a
    // fixed label list of size 6

    List<List<FixedList<label,6>>> globalData(Pstream::nProcs());

    forAll(cellsToRecv_, i)
    {
        if (cellsToRecv_[i].size() > 0)
        {
            FixedList<label,6> data(0);

            data[0] = neighbors_[i];
            data[1] = cellsToRecv_[i].size();
            data[2] = i;

            // Compute the boundary offset that the neighbor should have

            if (i > 0)
            {
                const labelVector bo =
                  - (neighborTs[i].T() & neighborOffsets_[i]);

                for (int j = 0; j < 3; j++)
                    data[j+3] = bo[j];
            }

            globalData[Pstream::myProcNo()].append(data);
        }
    }

    Pstream::gatherList(globalData, Pstream::msgType(), lvl.comms());
    Pstream::scatterList(globalData, Pstream::msgType(), lvl.comms());

    // Collect

    cellsToSend_.clear();
    cellsToSend_.setSize(lvl.boundaries().size()+1);

    sendTags_.clear();
    sendTags_.setSize(lvl.boundaries().size()+1, -1);

    forAll(globalData, proc)
    forAll(globalData[proc], i)
    {
        const FixedList<label,6>& data = globalData[proc][i];

        if (data[0] == Pstream::myProcNo())
        {
            // Requested boundary offset

            const labelVector bo(data[3], data[4], data[5]);

            // Find the matching neighbor

            bool found = false;

            forAll(neighbors_, j)
            if (neighbors_[j] == proc && bo == neighborOffsets_[j])
            {
                cellsToSend_[j].setSize(data[1]);
                sendTags_[j] = data[2];

                found = true;
                break;
            }

            if (!found)
                FatalErrorInFunction
                    << "Could not find matching neighbor for boundary offset "
                    << bo << endl << abort(FatalError);
        }
    }

    if (Pstream::parRun())
    {
        // Only exchange data when running in parallel

        forAll(cellsToSend_, i)
            if (cellsToSend_[i].size() > 0)
                UIPstream::read
                (
                    Pstream::commsTypes::nonBlocking,
                    neighbors_[i],
                    reinterpret_cast<char*>(cellsToSend_[i].begin()),
                    cellsToSend_[i].size()*sizeof(labelVector),
                    sendTags_[i],
                    lvl.comms()
                );

        forAll(cellsToRecv_, i)
            if (cellsToRecv_[i].size() > 0)
                UOPstream::write
                (
                    Pstream::commsTypes::nonBlocking,
                    neighbors_[i],
                    reinterpret_cast<char*>(cellsToRecv_[i].begin()),
                    cellsToRecv_[i].size()*sizeof(labelVector),
                    i,
                    lvl.comms()
                );

        Pstream::waitRequests();
    }
    else
    {
        // When running in serial, the only exchange can occur across periodic
        // boundaries. Find the matching send/recv cells manually based on the
        // neighbor offset.

        forAll(cellsToRecv_, i)
        if (cellsToRecv_[i].size() > 0)
        {
            bool found = false;

            forAll(cellsToSend_, j)
            if (cellsToSend_[j].size() > 0)
            if (neighborOffsets_[i] == -neighborOffsets_[j])
            {
                cellsToSend_[j] = cellsToRecv_[i];

                found = true;
                break;
            }

            if (!found)
                FatalErrorInFunction
                    << "Could not found periodic neighbor"
                    << abort(FatalError);
        }
    }

    // Displace and transform send cell indices. Not needed for local cells.

    for (label i = 1; i < cellsToSend_.size(); i++)
    {
        const labelTensor T = neighborTs[i];
        const faceLabel J = pTransform<faceLabel>(T,I);

        const labelVector LL = briscola::cmptMin(J.lower(), J.upper());
        const labelVector UU = briscola::cmptMax(J.lower(), J.upper());

        const labelVector bo = (T & neighborOffsets_[i]);

        forAll(cellsToSend_[i], j)
        {
            labelVector& index = cellsToSend_[i][j];

            for (int k = 0; k < 3; k++)
                if (bo[k] != 0)
                    index[k] += index[k] > 0 ? (LL[k]-1) : UU[k];

            // Transform back

            index = (T.T() & index);

            for (int k = 0; k < 3; k++)
                index[k] += index[k] < 0 ? UU[k] : 0;
        }
    }
}

template<class MeshType>
template<class Type>
List<Type> cellDataExchange<MeshType>::dataFunc
(
    const meshField<Type,MeshType>& field
) const
{
    List<Type> D(indices_.size(), Zero);

    // Prepare buffers

    List<List<Type>> recvData(cellsToRecv_.size());
    List<List<Type>> sendData(cellsToSend_.size());

    forAll(recvData, i)
        recvData[i].resize(cellsToRecv_[i].size());

    forAll(sendData, i)
        sendData[i].resize(cellsToSend_[i].size());

    forAll(sendData, i)
        forAll(sendData[i], j)
            sendData[i][j] = field(this->l_, this->d_, cellsToSend_[i][j]);

    // Only send/receive when running in parallel. When running in serial, copy
    // manually

    if (Pstream::parRun())
    {
        // Receive

        forAll(recvData, i)
            if (recvData[i].size() > 0)
                UIPstream::read
                (
                    Pstream::commsTypes::nonBlocking,
                    neighbors_[i],
                    reinterpret_cast<char*>(recvData[i].begin()),
                    recvData[i].byteSize(),
                    i,
                    field[this->l_].lvl().comms()
                );

        // Send

        forAll(sendData, i)
            if (sendData[i].size() > 0)
                UOPstream::write
                (
                    Pstream::commsTypes::nonBlocking,
                    neighbors_[i],
                    reinterpret_cast<char*>(sendData[i].begin()),
                    sendData[i].byteSize(),
                    sendTags_[i],
                    field[this->l_].lvl().comms()
                );

        Pstream::waitRequests();
    }
    else
    {
        forAll(recvData, i)
            if (recvData[i].size() > 0)
                forAll(sendData, j)
                    if (sendData[j].size() > 0)
                        if (neighborOffsets_[i] == -neighborOffsets_[j])
                            recvData[i] = sendData[j];
    }

    // Copy back from buffer to list

    forAll(D, i)
        D[i] = recvData[map_[i].first()][map_[i].second()];

    return D;
}

// Instantiate class and data functions

template class cellDataExchange<colocated>;
template class cellDataExchange<staggered>;

#define INSTANTIATE(TYPE,MESHTYPE)                                             \
template List<TYPE> cellDataExchange<MESHTYPE>::dataFunc                       \
(                                                                              \
    const meshField<TYPE,MESHTYPE>&                                            \
) const;

INSTANTIATE(scalar,colocated)
INSTANTIATE(vector,colocated)
INSTANTIATE(tensor,colocated)
INSTANTIATE(sphericalTensor,colocated)
INSTANTIATE(symmTensor,colocated)
INSTANTIATE(diagTensor,colocated)
INSTANTIATE(faceScalar,colocated)
INSTANTIATE(faceVector,colocated)

INSTANTIATE(scalar,staggered)
INSTANTIATE(vector,staggered)
INSTANTIATE(tensor,staggered)
INSTANTIATE(sphericalTensor,staggered)
INSTANTIATE(symmTensor,staggered)
INSTANTIATE(diagTensor,staggered)
INSTANTIATE(faceScalar,staggered)
INSTANTIATE(faceVector,staggered)

#undef INSTANTIATE

}

}

}
