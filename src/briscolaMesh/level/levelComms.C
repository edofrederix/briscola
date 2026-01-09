#include "levelComms.H"

#include "level.H"
#include "mesh.H"

#include "PstreamGlobals.H"
#include "parallelBoundary.H"

namespace Foam
{

namespace briscola
{

MPI_Comm& levelComms::communicator()
{
    return PstreamGlobals::MPICommunicators_[comm_];
}

levelComms::levelComms(const level& lvl)
:
    lvl_(lvl),
    comm_(-1)
{
    if (Pstream::parRun())
    {
        DynamicList<label> neighbors, weights;
        DynamicList<DynamicList<label>> sendTags;

        forAll(lvl_.boundaries(), i)
        if (lvl_.boundaries()[i].castable<parallelBoundary>())
        {
            const parallelBoundary& b =
                lvl_.boundaries()[i].cast<parallelBoundary>();

            const labelVector bo(b.offset());
            const labelVector N(lvl_.N());

            const label proc = b.neighborProcNum();
            const label tag = b.tag();

            // Do not handle periodic boundaries to self

            if (proc == Pstream::myProcNo())
                continue;

            label weight = 1;

            for (int d = 0; d < 3; d++)
                if (bo[d] == 0)
                    weight *= N[d];

            const label index = findIndex(neighbors, proc);

            if (index != -1)
            {
                weights[index] += weight;
                sendTags[index].append(tag);
            }
            else
            {
                neighbors.append(proc);
                weights.append(weight);
                sendTags.append(DynamicList<label>(1,tag));
            }
        }

        // Create a new communicator using OpenFOAM's methodology so that it is
        // registered

        labelList validProcs(Pstream::nProcs());
        forAll(validProcs, proc)
            validProcs[proc] = proc;

        comm_ = Pstream::allocateCommunicator(UPstream::worldComm, validProcs);

        MPI_Comm& comm = communicator();

        // Update the communicator to a neighbor graph communicator, which can
        // be used by MPI neighbor exchange functions

        MPI_Dist_graph_create_adjacent
        (
            PstreamGlobals::MPI_COMM_FOAM,
            neighbors.size(),
            neighbors.begin(),
            weights.begin(),
            neighbors.size(),
            neighbors.begin(),
            weights.begin(),
            MPI_INFO_NULL,
            true,
            &comm
        );

        // Communicator checks

        label indeg, outdeg, weighted;
        MPI_Dist_graph_neighbors_count
        (
            comm,
            &indeg,
            &outdeg,
            &weighted
        );

        // In- and out-degree must match the number of neighbors

        if (neighbors.size() != indeg)
            FatalError
                << "Error setting up graph communicator"
                << abort(FatalError);

        if (neighbors.size() != outdeg)
            FatalError
                << "Error setting up graph communicator"
                << abort(FatalError);

        labelList srcweights(indeg);
        labelList dstweights(outdeg);

        recvNeighbors_.resize(indeg);
        sendNeighbors_.resize(outdeg);

        MPI_Dist_graph_neighbors
        (
            comm,
            indeg,
            recvNeighbors_.begin(),
            srcweights.begin(),
            outdeg,
            sendNeighbors_.begin(),
            dstweights.begin()
        );

        // Communicate tags

        DynamicList<label> sendBuffer;
        forAll(sendTags, i)
            forAll(sendTags[i], j)
                sendBuffer.append(sendTags[i][j]);

        labelList recvBuffer(sendBuffer.size());

        labelList sendCount(outdeg);
        labelList recvCount(indeg);

        forAll(sendCount, i)
            sendCount[i] = sendTags[i].size();

        forAll(recvCount, i)
            recvCount[findIndex(recvNeighbors_,sendNeighbors_[i])] =
                sendTags[i].size();

        labelList sendDisp(outdeg, 0);
        labelList recvDisp(indeg, 0);

        for(label i = 1; i < outdeg; i++)
            sendDisp[i] = sendDisp[i-1] + sendCount[i-1];

        for(label i = 1; i < indeg; i++)
            recvDisp[i] = recvDisp[i-1] + recvCount[i-1];

        MPI_Neighbor_alltoallv
        (
            sendBuffer.begin(),
            sendCount.begin(),
            sendDisp.begin(),
            MPI_INT,
            recvBuffer.begin(),
            recvCount.begin(),
            recvDisp.begin(),
            MPI_INT,
            comm
        );

        // Collect receive tags

        List<DynamicList<label>> recvTags(indeg);

        forAll(recvTags, i)
            for (label j = 0; j < recvCount[i]; j++)
                recvTags[i].append(recvBuffer[recvDisp[i] + j]);

        // Set send map

        sendMap_.resize(outdeg);

        forAll(sendNeighbors_, i)
        {
            const label proc = sendNeighbors_[i];

            sendMap_[i].resize(sendTags[i].size());
            sendMap_[i] = -1;

            forAll(sendTags[i], j)
            {
                const label tag = sendTags[i][j];

                // Find corresponding boundary for given proc and tag

                bool found = false;

                forAll(lvl_.boundaries(), bi)
                if (lvl_.boundaries()[bi].castable<parallelBoundary>())
                {
                    const parallelBoundary& b =
                        lvl_.boundaries()[bi].cast<parallelBoundary>();

                    if (b.neighborProcNum() == proc && b.tag() == tag)
                    {
                        sendMap_[i][j] = bi;
                        found = true;
                        break;
                    }
                }

                if (!found)
                    FatalErrorInFunction
                        << "Could not find matching proc/tag pair"
                        << abort(FatalError);
            }
        }

        // Set receive map

        recvMap_.resize(indeg);

        forAll(recvNeighbors_, i)
        {
            const label proc = recvNeighbors_[i];

            recvMap_[i].resize(recvTags[i].size());

            forAll(recvTags[i], j)
            {
                const label tag = recvTags[i][j];

                // Find corresponding boundary for given proc and tag

                bool found = false;

                forAll(lvl_.boundaries(), bi)
                if (lvl_.boundaries()[bi].castable<parallelBoundary>())
                {
                    const parallelBoundary& b =
                        lvl_.boundaries()[bi].cast<parallelBoundary>();

                    if (b.neighborProcNum() == proc && b.tag() == tag)
                    {
                        recvMap_[i][j] = bi;
                        found = true;
                        break;
                    }
                }

                if (!found)
                    FatalErrorInFunction
                        << "Could not find matching proc/tag pair"
                        << abort(FatalError);
            }
        }
    }
}

levelComms::levelComms(const levelComms& c)
:
    lvl_(c.lvl_),
    comm_(c.comm_),
    recvNeighbors_(c.recvNeighbors_),
    sendNeighbors_(c.sendNeighbors_),
    recvMap_(c.recvMap_),
    sendMap_(c.sendMap_)
{}

levelComms::levelComms(const levelComms& c, const level& lvl)
:
    lvl_(lvl),
    comm_(c.comm_),
    recvNeighbors_(c.recvNeighbors_),
    sendNeighbors_(c.sendNeighbors_),
    recvMap_(c.recvMap_),
    sendMap_(c.sendMap_)
{}

levelComms::~levelComms()
{}

const MPI_Comm& levelComms::communicator() const
{
    return PstreamGlobals::MPICommunicators_[comm_];
}

}

}
