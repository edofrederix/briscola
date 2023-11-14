#include "cellDataExchange.H"
#include "boundaryPartPatch.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTemplateTypeNameAndDebug(cellDataExchange<colocated>, 0);
defineTemplateTypeNameAndDebug(cellDataExchange<staggered>, 0);

template<class MeshType>
void cellDataExchange<MeshType>::setupComms
(
    const labelList& myProcsToTalkTo,
    const labelList& recvCounts
)
{
    recvProcNums_ = myProcsToTalkTo;
    recvCounts_ = recvCounts;

    List<List<Tuple2<label,label>>> globalData(Pstream::nProcs());

    globalData[Pstream::myProcNo()].setSize(myProcsToTalkTo.size());

    forAll(myProcsToTalkTo, i)
        globalData[Pstream::myProcNo()][i] =
            Tuple2<label,label>(myProcsToTalkTo[i], recvCounts[i]);

    Pstream::gatherList(globalData);
    Pstream::scatterList(globalData);

    // Collect

    sendProcNums_.setSize(0);
    sendCounts_.setSize(0);

    forAll(globalData, i)
    if (i != Pstream::myProcNo())
    {
        forAll(globalData[i], j)
        {
            if (globalData[i][j].first() == Pstream::myProcNo())
            {
                sendProcNums_.append(i);
                sendCounts_.append(globalData[i][j].second());
            }
        }
    }

    // Set patch pointers if they exist

    recvPatchPtrs_.setSize(recvProcNums_.size());
    sendPatchPtrs_.setSize(sendProcNums_.size());

    recvPatchPtrs_ = nullptr;
    sendPatchPtrs_ = nullptr;

    const mesh& msh = this->fvMsh_.msh();

    forAll(msh.partPatches(), i)
    {
        const partPatch& patch = msh.partPatches()[i];

        if (patch.typeNum() != boundaryPartPatch::typeNumber)
        {
            const label procNum = patch.neighborProcNum();

            const label ir = findIndex(recvProcNums_, procNum);
            const label is = findIndex(sendProcNums_, procNum);

            if (ir > -1)
                sendPatchPtrs_[ir] = &patch;

            if (is > -1)
                recvPatchPtrs_[is] = &patch;
        }
    }
}

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
    sendCells_(e.sendCells_),
    recvCells_(e.recvCells_),
    sendProcNums_(e.sendProcNums_),
    recvProcNums_(e.recvProcNums_),
    sendCounts_(e.sendCounts_),
    recvCounts_(e.recvCounts_),
    sendPatchPtrs_(e.sendPatchPtrs_),
    recvPatchPtrs_(e.recvPatchPtrs_)
{}

template<class MeshType>
cellDataExchange<MeshType>::~cellDataExchange()
{}

template<class MeshType>
void cellDataExchange<MeshType>::init(const List<labelVector>& indices)
{
    indices_ = indices;

    map_.clear();
    map_.setSize(indices.size());

    const faceLabel I = this->fvMsh_.template I<MeshType>(this->l_,this->d_);

    const labelVector U = I.upper();
    const labelVector L = I.lower();

    const mesh& msh = this->fvMsh_.msh();

    // Collect the required cells and store per neighbor processor

    recvCells_.clear();
    recvCells_.setSize(msh.partPatches().size());

    forAll(indices, i)
    {
        labelVector index = indices[i];

        const labelVector closest = labelVector
        (
            Foam::max(Foam::min(index.x(), U.x()-1), L.x()),
            Foam::max(Foam::min(index.y(), U.y()-1), L.y()),
            Foam::max(Foam::min(index.z(), U.z()-1), L.z())
        );

        const labelVector offset = index - closest;

        if (offset == zeroXYZ)
            FatalErrorInFunction
                << "Requesting data exchange for a local cell" << endl
                << abort(FatalError);

        const labelVector bo =
            briscola::cmptMin(briscola::cmptMax(offset, -unitXYZ), unitXYZ);

        const label degree = cmptSum(cmptMag(bo));

        if (!this->fvMsh_.structured() && degree > 1)
            FatalErrorInFunction
                << "Can only obtain cell data across faces "
                << "on unstructured meshes" << endl
                << abort(FatalError);

        // Get the patch to the other processor

        label patchNum = -1;

        forAll(msh.partPatches(), patchi)
        if (msh.partPatches()[patchi].boundaryOffset() == bo)
        {
            patchNum = patchi;
            break;
        }

        if (patchNum == -1)
            FatalErrorInFunction
                << "Could not find patch corresponding to boundary offset "
                << bo << endl << abort(FatalError);

        const partPatch& patch = msh.partPatches()[patchNum];

        if (patch.typeNum() == boundaryPartPatch::typeNumber)
            FatalErrorInFunction
                << "Cannot exchange cell data across a boundary part patch"
                << endl << abort(FatalError);

        // Store trimmed index

        for (int j = 0; j < 3; j++)
            if (bo[j] != 0)
                index[j] = bo[j] > 0 ? index[j]-(U[j]-1) : index[j]-L[j];

        recvCells_[patchNum].append(index);

        map_[i] =
            Tuple2<label,label>
            (
                patchNum,
                recvCells_[patchNum].size()-1
            );
    }

    // Exchange processor numbers to talk to

    labelList myProcsToTalkTo;
    labelList cellCount;

    forAll(recvCells_, i)
    if (recvCells_[i].size() > 0)
    {
        myProcsToTalkTo.append(msh.partPatches()[i].neighborProcNum());
        cellCount.append(recvCells_[i].size());
    }

    setupComms(myProcsToTalkTo, cellCount);

    // Trim recv data

    List<List<labelVector>> recvCellsTemp;

    forAll(recvCells_, i)
        if (recvCells_[i].size() > 0)
            recvCellsTemp.append(recvCells_[i]);

    recvCells_ = recvCellsTemp;

    // Send indices to neighbors and receive from neighbors

    sendCells_.clear();
    sendCells_.setSize(sendProcNums_.size());

    forAll(sendProcNums_, i)
    {
        sendCells_[i].setSize(sendCounts_[i]);

        UIPstream::read
        (
            Pstream::commsTypes::nonBlocking,
            sendProcNums_[i],
            reinterpret_cast<char*>(sendCells_[i].begin()),
            sendCells_[i].byteSize()
        );
    }

    forAll(recvProcNums_, i)
    {
        UOPstream::write
        (
            Pstream::commsTypes::nonBlocking,
            recvProcNums_[i],
            reinterpret_cast<char*>(recvCells_[i].begin()),
            recvCells_[i].byteSize()
        );
    }

    UPstream::waitRequests();

    // Transform send cell indices

    forAll(sendProcNums_, i)
    {
        const partPatch& patch = *sendPatchPtrs_[i];
        const labelTensor& T = patch.T();

        faceLabel J = pTransform<faceLabel>(T,I);

        const labelVector L = briscola::cmptMin(J.lower(), J.upper());
        const labelVector U = briscola::cmptMax(J.lower(), J.upper());

        const labelVector bo = (T & patch.boundaryOffset());

        forAll(sendCells_[i], j)
        {
            labelVector& index = sendCells_[i][j];

            for (int k = 0; k < 3; k++)
                if (bo[k] != 0)
                    index[k] = index[k] > 0 ? (L[k]-1)+index[k] : U[k]+index[k];

            // Transform back

            index = (T.T() & index);

            for (int k = 0; k < 3; k++)
                index[k] = index[k] < 0 ? U[k]+index[k] : index[k];
        }
    }
}

template<class MeshType>
template<class Type>
List<Type> cellDataExchange<MeshType>::dataFunc
(
    const meshField<Type,MeshType>& field
)
{
    List<List<Type>> data(recvCells_.size());

    forAll(recvCells_, i)
    {
        data[i].setSize(recvCounts_[i]);

        UIPstream::read
        (
            Pstream::commsTypes::nonBlocking,
            recvProcNums_[i],
            reinterpret_cast<char*>(data[i].begin()),
            data[i].byteSize()
        );
    }

    forAll(sendCells_, i)
    {
        List<Type> data(sendCells_[i].size());

        forAll(data, j)
            data[j] = field(this->l_,this->d_,sendCells_[i][j]);

        UOPstream::write
        (
            Pstream::commsTypes::nonBlocking,
            sendProcNums_[i],
            reinterpret_cast<char*>(data.begin()),
            data.byteSize()
        );
    }

    UPstream::waitRequests();

    List<Type> D(indices_.size());

    forAll(D, i)
        D[i] = data[map_[i].first()][map_[i].second()];

    return D;
}

// Instantiate class and data functions

template class cellDataExchange<colocated>;
template class cellDataExchange<staggered>;

#define INSTANTIATE(TYPE,MESHTYPE)                                             \
template List<TYPE> cellDataExchange<MESHTYPE>::dataFunc                       \
(                                                                              \
    const meshField<TYPE,MESHTYPE>&                                            \
);

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
