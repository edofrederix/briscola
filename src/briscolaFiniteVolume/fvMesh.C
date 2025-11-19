#include "fvMesh.H"

#include "colocatedFields.H"
#include "staggeredFields.H"
#include "immersedBoundary.H"

#include "PstreamGlobals.H"
#include "parallelBoundary.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(fvMesh, 0);

template<class MeshType>
void fvMesh::setInternalCells()
{
    faceLabel slave;
    for (int i = 0; i < 6; i++)
        slave[i] = mshPtr_->b(faceOffsets[i]).slave();

    forAll(*this, l)
    {
        const part& p = this->operator[](l);

        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            faceLabel& I = this->I<MeshType>(l,d);

            const labelVector& padding = MeshType::padding[d];

            I = faceLabel(zeroXYZ, p.N()+padding);

            for (int i = 0; i < 6; i++)
                I[i] += (padding[i/2] && slave[i]) ? 1 - 2*(i%2) : 0;
        }
    }
}

template<>
faceLabel& fvMesh::I<colocated>(const label l, const label d)
{
    return Ic_[l][d];
}

template<>
faceLabel& fvMesh::I<staggered>(const label l, const label d)
{
    return Is_[l][d];
}

void fvMesh::setDistributedCommGraph()
{
    if (Pstream::parRun() && MPI_COMM_FOAM_OLD == MPI_COMM_NULL)
    {
        DynamicList<label> neighbors, weights;
        DynamicList<DynamicList<label>> sendTags;

        forAll(mshPtr_->boundaries(), i)
        if (mshPtr_->boundaries()[i].castable<parallelBoundary>())
        {
            const parallelBoundary& b =
                mshPtr_->boundaries()[i].cast<parallelBoundary>();

            const labelVector bo(b.offset());
            const labelVector N(mshPtr_->operator[](0).N());

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

        // Overwrite the default OpenFOAM communicator with a neighbor graph
        // communicator, which can be used by MPI neighbor exchange functions.
        // This new and optimized communicator is then automatically used by all
        // operations including reductions and gFunctions.

        MPI_COMM_FOAM_OLD = PstreamGlobals::MPI_COMM_FOAM;

        MPI_Dist_graph_create_adjacent
        (
            MPI_COMM_FOAM_OLD,
            neighbors.size(),
            neighbors.begin(),
            weights.begin(),
            neighbors.size(),
            neighbors.begin(),
            weights.begin(),
            MPI_INFO_NULL,
            true,
            &PstreamGlobals::MPI_COMM_FOAM
        );

        // Communicator checks

        label indeg, outdeg, weighted;
        MPI_Dist_graph_neighbors_count
        (
            PstreamGlobals::MPI_COMM_FOAM,
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

        MPI_neighbors_recv.resize(indeg);
        MPI_neighbors_send.resize(outdeg);

        MPI_Dist_graph_neighbors
        (
            PstreamGlobals::MPI_COMM_FOAM,
            indeg,
            MPI_neighbors_recv.begin(),
            srcweights.begin(),
            outdeg,
            MPI_neighbors_send.begin(),
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
            recvCount[findIndex(MPI_neighbors_recv,MPI_neighbors_send[i])] =
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
            PstreamGlobals::MPI_COMM_FOAM
        );

        // Collect receive tags

        List<DynamicList<label>> recvTags(indeg);

        forAll(recvTags, i)
            for (label j = 0; j < recvCount[i]; j++)
                recvTags[i].append(recvBuffer[recvDisp[i] + j]);

        // Set send map

        MPI_neighbors_send_map.resize(outdeg);

        forAll(MPI_neighbors_send, i)
        {
            const label proc = MPI_neighbors_send[i];

            MPI_neighbors_send_map[i].resize(sendTags[i].size());
            MPI_neighbors_send_map[i] = -1;

            forAll(sendTags[i], j)
            {
                const label tag = sendTags[i][j];

                // Find corresponding boundary for given proc and tag

                bool found = false;

                forAll(mshPtr_->boundaries(), bi)
                if (mshPtr_->boundaries()[bi].castable<parallelBoundary>())
                {
                    const parallelBoundary& b =
                        mshPtr_->boundaries()[bi].cast<parallelBoundary>();

                    if (b.neighborProcNum() == proc && b.tag() == tag)
                    {
                        MPI_neighbors_send_map[i][j] = bi;
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

        MPI_neighbors_recv_map.resize(indeg);

        forAll(MPI_neighbors_recv, i)
        {
            const label proc = MPI_neighbors_recv[i];

            MPI_neighbors_recv_map[i].resize(recvTags[i].size());

            forAll(recvTags[i], j)
            {
                const label tag = recvTags[i][j];

                // Find corresponding boundary for given proc and tag

                bool found = false;

                forAll(mshPtr_->boundaries(), bi)
                if (mshPtr_->boundaries()[bi].castable<parallelBoundary>())
                {
                    const parallelBoundary& b =
                        mshPtr_->boundaries()[bi].cast<parallelBoundary>();

                    if (b.neighborProcNum() == proc && b.tag() == tag)
                    {
                        MPI_neighbors_recv_map[i][j] = bi;
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

fvMesh::fvMesh(const IOdictionary& dict, const Time& time)
:
    regIOobject(dict, true),
    mshPtr_(mesh::New(dict)),
    schemeDict_
    (
        IOobject
        (
            "briscolaSchemeDict",
            time.system(),
            time,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    solverDict_
    (
        IOobject
        (
            "briscolaSolverDict",
            time.system(),
            time,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    Ic_(mshPtr_->size(), List<faceLabel>(colocated::numberOfDirections)),
    Is_(mshPtr_->size(), List<faceLabel>(staggered::numberOfDirections)),
    colocatedMetrics_(),
    staggeredMetrics_()
{
    setInternalCells<colocated>();
    setInternalCells<staggered>();

    setDistributedCommGraph();

    colocatedMetrics_ = new fvMeshMetrics<colocated>(*this);

    // Only generate staggered metrics when the brick topology is structured

    if (mshPtr_->structured())
        staggeredMetrics_ = new fvMeshMetrics<staggered>(*this);
}

fvMesh::fvMesh(const fvMesh& fvMsh)
:
    regIOobject(fvMsh, true),
    mshPtr_(fvMsh.mshPtr_, false),
    schemeDict_(fvMsh.schemeDict_),
    solverDict_(fvMsh.solverDict_),
    Ic_(fvMsh.Ic_),
    Is_(fvMsh.Is_),
    colocatedMetrics_(),
    staggeredMetrics_()
{
    colocatedMetrics_ = new fvMeshMetrics<colocated>(*this);

    // Only generate staggered metrics when the brick topology is structured

    if (mshPtr_->structured())
        staggeredMetrics_ = new fvMeshMetrics<staggered>(*this);
}

fvMesh::~fvMesh()
{
    if (MPI_COMM_FOAM_OLD != MPI_COMM_NULL)
        MPI_Comm_free(&MPI_COMM_FOAM_OLD);
}

template<>
const faceLabel& fvMesh::I<colocated>(const label l, const label d) const
{
    return Ic_[l][d];
}

template<>
const faceLabel& fvMesh::I<staggered>(const label l, const label d) const
{
    return Is_[l][d];
}

template<>
const fvMeshMetrics<colocated>& fvMesh::metrics<colocated>() const
{
    return colocatedMetrics_();
}

template<>
const fvMeshMetrics<staggered>& fvMesh::metrics<staggered>() const
{
    #ifdef FULLDEBUG
    if (!mshPtr_->structured())
    {
        FatalErrorInFunction
            << "Staggered metrics are not generated on unstructured meshes."
            << endl << abort(FatalError);
    }
    #endif

    return staggeredMetrics_();
}

template<>
labelVector fvMesh::findCell<colocated>
(
    const vector& p,
    const label l,
    const label
) const
{
    return mshPtr_->findCell(p,l);
}

template<>
labelVector fvMesh::findCell<staggered>
(
    const vector& p,
    const label l,
    const label d
) const
{
    labelVector colo = mshPtr_->findCell(p,l);

    const meshDirection<vertexVector,staggered>& v =
        this->metrics<staggered>().vertexCenters()[l][d];

    labelVector vU(v.I().upper());

    if (colo != -unitXYZ)
    {
        // Search staggered cells that may overlap with the colocated one

        labelVector L(colo);
        labelVector U(colo + staggered::padding[d] + unitXYZ);

        if (this->rectilinear() != unitXYZ)
        {
            // Add surrounding cells in non-shifted direction

            for (int dir = 0; dir < 3; dir++)
            if (dir != d)
            {
                L -= staggered::padding[dir];
                U += staggered::padding[dir];
            }

            // Trim

            for (int dir = 0; dir < 3; dir++)
            if (dir != d)
            {
                L[dir] = Foam::max(L[dir], 0);
                U[dir] = Foam::min(U[dir], vU[dir]);
            }
        }

        for (int i = L.x(); i < U.x(); i++)
            for (int j = L.y(); j < U.y(); j++)
                for (int k = L.z(); k < U.z(); k++)
                    if (interpolationWeights(p,v(i,j,k),true) != -vector::one)
                        return labelVector(i,j,k);

        // When the point is near a shifted boundary, it may be found on the
        // colocated mesh but not on the staggered mesh.

        return -unitXYZ;
    }
    else
    {
        // Near shifted boundaries, the point may be outside the colocated mesh
        // but still on the staggered mesh. This happens in two situations:
        //
        // 1) The point is outside the physical domain -> return -unitXYZ
        // 2) The point is on a parallel/periodic master/slave slab of cells.
        //    Data on such cells is redundantly stored across two processors, so
        //    we let the processor which has the point on its colocated mesh
        //    handle the situation -> return -unitXYZ

        return -unitXYZ;
    }
}

template<class MeshType>
List<labelVector> fvMesh::findCells
(
    const vectorList& points,
    const label l,
    const label d
) const
{
    List<labelVector> res(points.size());

    forAll(points, i)
        res[i] = this->findCell<MeshType>(points[i], l, d);

    return res;
}

// Instantiate

template List<labelVector>
fvMesh::findCells<colocated>(const vectorList&, const label, const label) const;
template List<labelVector>
fvMesh::findCells<staggered>(const vectorList&, const label, const label) const;

template const FastPtrList<immersedBoundary<staggered>>&
fvMesh::immersedBoundaries<staggered>() const;
template const FastPtrList<immersedBoundary<colocated>>&
fvMesh::immersedBoundaries<colocated>() const;

MPI_Comm fvMesh::MPI_COMM_FOAM_OLD = MPI_COMM_NULL;

labelList fvMesh::MPI_neighbors_recv;
labelList fvMesh::MPI_neighbors_send;

List<labelList> fvMesh::MPI_neighbors_recv_map;
List<labelList> fvMesh::MPI_neighbors_send_map;

}

}

}
