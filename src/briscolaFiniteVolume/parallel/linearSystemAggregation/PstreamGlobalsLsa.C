#include "PstreamGlobalsLsa.H"

namespace Foam
{

namespace PstreamGlobals
{

DynamicList<Tuple2<label,label>> lsaCommsData_;
DynamicList<Tuple2<label,label>> lsaMasterCommsData_;

label lsaSetComms
(
    const label nParts,
    const label color,
    const label key,
    const bool masters
)
{
    if (!Pstream::parRun())
        return -1;

    DynamicList<Tuple2<label,label>>& lsaCommsData =
        masters ? lsaMasterCommsData_ : lsaCommsData_;

    // Check if a suitable communicator already exists

    forAll(lsaCommsData, i)
        if (lsaCommsData[i].first() == nParts)
            return lsaCommsData[i].second();

    // Create new one otherwise

    MPICommunicators_.append(MPI_Comm());

    const label commNum = MPICommunicators_.size() - 1;

    MPI_Comm_split
    (
        MPI_COMM_FOAM,
        color > -1 ? color: MPI_UNDEFINED,
        key,
        &MPICommunicators_[commNum]
    );

    lsaCommsData.append(Tuple2<label,label>(nParts, commNum));

    return commNum;
}

MPI_Comm& lsaGetComm(const label commNum_)
{
    return
        commNum_ < 0
      ? PstreamGlobals::MPI_COMM_FOAM
      : MPICommunicators_[commNum_];
}

}

}
