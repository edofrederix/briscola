#include "PstreamSplit.H"

namespace Foam
{

namespace PstreamGlobals
{

DynamicList<MPI_Comm> splitComms_;
DynamicList<List<bool>> splitCommMembers_;

label split(const bool member)
{
    if (!Pstream::parRun())
        return -1;

    List<bool> members(Pstream::nProcs());
    members[Pstream::myProcNo()] = member;

    Pstream::gatherList(members);
    Pstream::scatterList(members);

    // Check if a suitable communicator already exists

    forAll(splitCommMembers_, comm)
        if (splitCommMembers_[comm] == members)
            return comm;

    // Create new one otherwise

    splitComms_.append(MPI_Comm());
    splitCommMembers_.append(members);

    const label comm = splitComms_.size() - 1;

    MPI_Comm_split
    (
        PstreamGlobals::MPI_COMM_FOAM,
        member ? 1 : MPI_UNDEFINED,
        Pstream::myProcNo(),
        &splitComms_[comm]
    );

    return comm;
}

MPI_Comm& getSplitComm(const label comm)
{
    return splitComms_[comm];
}

}

}
