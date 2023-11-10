#include "dataExchange.H"
#include "boundaryPartPatch.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTemplateTypeNameAndDebug(dataExchange<colocated>, 0);
defineTemplateTypeNameAndDebug(dataExchange<staggered>, 0);

template<class MeshType>
void dataExchange<MeshType>::setupComms
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

    const mesh& msh = fvMsh_.msh();

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
dataExchange<MeshType>::dataExchange
(
    const fvMesh& fvMsh,
    const label d
)
:
    fvMsh_(fvMsh),
    d_(d)
{}

template<class MeshType>
dataExchange<MeshType>::dataExchange
(
    const dataExchange& e
)
:
    fvMsh_(e.fvMsh_),
    d_(e.d_),
    sendProcNums_(e.sendProcNums_),
    recvProcNums_(e.recvProcNums_),
    sendCounts_(e.sendCounts_),
    recvCounts_(e.recvCounts_),
    sendPatchPtrs_(e.sendPatchPtrs_),
    recvPatchPtrs_(e.recvPatchPtrs_)
{}

template<class MeshType>
dataExchange<MeshType>::~dataExchange()
{}

template class dataExchange<colocated>;
template class dataExchange<staggered>;

}

}

}
