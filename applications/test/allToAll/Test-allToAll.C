#include "arguments.H"

#include "vector.H"
#include "block.H"

#include <chrono>

using namespace Foam;
using namespace briscola;

// Return the processor number given an index and decomposition

label procNumFromIndex
(
    const labelVector ijk,
    const labelVector D
)
{
    return ijk.x()*D.y()*D.z() + ijk.y()*D.z() + ijk.z();
}

// Return the index given a processor number and decomposition

labelVector indexFromProcNum
(
    const label num,
    const labelVector D
)
{
    label i = num/(D.y()*D.z());
    label j = (num-i*D.y()*D.z())/D.z();
    label k = (num-i*D.y()*D.z()-j*D.z());

    return labelVector(i,j,k);
}

// Function to unpack recv buffer

void unpack
(
    const vectorBlock& buffer,
    const List<labelVector>& recvStart,
    const List<labelVector>& recvSize,
    vectorBlock& output,
    List<labelVector>& startIndex
)
{
    label cursor = 0;

    labelVector N(output.shape());

    for (label proc = 0; proc < Pstream::nProcs(); proc++)
    {
        labelVector start
        (
            recvStart[proc].x() - startIndex[Pstream::myProcNo()].x(),
            recvStart[proc].y() - startIndex[Pstream::myProcNo()].y(),
            recvStart[proc].z() - startIndex[Pstream::myProcNo()].z()
        );

        labelVector size = recvSize[proc];

        for (int i = start.x(); i < start.x() + size.x(); i++)
        for (int j = start.y(); j < start.y() + size.y(); j++)
        for (int k = start.z(); k < start.z() + size.z(); k++)
        {
            output(i,j,k) = buffer(cursor++);
        }
    }
}

// Function to test if the data is correct

bool check(const vectorBlock& data, const labelVector offset)
{
    forAllBlock(data, i, j, k)
    {
        labelVector test = offset + labelVector(i,j,k);

        if (labelVector(data(i,j,k)) != test)
            FatalError
                << "Test failed. Expected " << test
                << " but received " << data(i,j,k) << endl
                << "At index " << labelVector(i,j,k) << endl
                << abort(FatalError);
    }

    return true;
}

int main(int argc, char *argv[])
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    arguments::addBoolOption("parallel", "run in parallel");

    arguments::validArgs.append("mesh size");
    arguments::validArgs.append("initial decomposition");
    arguments::validArgs.append("target decomposition");
    arguments::validArgs.append("number of send repetitions");

    arguments args(argc, argv);

    const labelVector N(args.argRead<labelVector>(1));
    const labelVector I(args.argRead<labelVector>(2));
    const labelVector T(args.argRead<labelVector>(3));
    const label repetitions(args.argRead<label>(4));

    // Some checks

    if (!Pstream::parRun())
    {
        FatalError << "Must be run in parallel" << endl;
        FatalError.exit();
    }

    if (cmptProduct(I) != Pstream::nProcs())
    {
        FatalError
            << "Initial decomposition does not match the number of processors"
            << endl;
        FatalError.exit();
    }

    if (cmptProduct(T) != Pstream::nProcs())
    {
        FatalError
            << "Target decomposition does not match the number of processors"
            << endl;
        FatalError.exit();
    }

    for (int dir = 0; dir < 3; dir++)
    if (N[dir] < I[dir])
    {
        FatalError
            << "Mesh not decomposable by initial decomposition vector"
            << endl;
        FatalError.exit();
    }

    for (int dir = 0; dir < 3; dir++)
    if (N[dir] < T[dir])
    {
        FatalError
            << "Mesh not decomposable by target decomposition vector"
            << endl;
        FatalError.exit();
    }

    // Print initial and target processor topologies

    Info<< "Initial topology = " << nl;

    for (label proc = 0; proc < Pstream::nProcs(); proc++)
        Info<< "    Proc " << proc << " at "
            << indexFromProcNum(proc, I) << endl;

    Info<< nl;
    Info<< "Target topology = " << endl;

    for (label proc = 0; proc < Pstream::nProcs(); proc++)
        Info<< "    Proc " << proc << " at "
            << indexFromProcNum(proc, T) << endl;

    Info<< endl;

    // Compute processor overlap from initial to target, defined by the global
    // start cell index and block size

    // Number of cells in initial and target decompositions
    List<labelVector> Ni(Pstream::nProcs(),Zero);
    List<labelVector> Nt(Pstream::nProcs(),Zero);

    // Start cell indices of each processor in initial and target decompositions
    List<labelVector> si(Pstream::nProcs(),Zero);
    List<labelVector> st(Pstream::nProcs(),Zero);
    for ( int proc = 0; proc < Pstream::nProcs(); proc++ )
    {
        Ni[proc] = vector
        (
            floor(N.x()/I.x()),
            floor(N.y()/I.y()),
            floor(N.z()/I.z())
        );

        // If the data distribution is uneven, remainder > 0
        labelVector remainder = N - cmptMultiply(Ni[proc], I);
        if (indexFromProcNum(proc, I).x() < remainder.x())
        {
            Ni[proc].x() += 1;
        }
        if (indexFromProcNum(proc, I).y() < remainder.y())
        {
            Ni[proc].y() += 1;
        }
        if (indexFromProcNum(proc, I).z() < remainder.z())
        {
            Ni[proc].z() += 1;
        }

        Nt[proc] = vector
        (
            floor(N.x()/T.x()),
            floor(N.y()/T.y()),
            floor(N.z()/T.z())
        );

        remainder = N - cmptMultiply(Nt[proc], T);
        if (indexFromProcNum(proc, T).x() < remainder.x())
        {
            Nt[proc].x() += 1;
        }
        if (indexFromProcNum(proc, T).y() < remainder.y())
        {
            Nt[proc].y() += 1;
        }
        if (indexFromProcNum(proc, T).z() < remainder.z())
        {
            Nt[proc].z() += 1;
        }

        // Set starting indices of each processor, initial decomp
        si[proc] = labelVector(0,0,0);
        for ( int x = 0; x < indexFromProcNum(proc, I).x(); x++ )
        {
            si[proc].x() += Ni[procNumFromIndex(vector(x,0,0), I)].x();
        }
        for ( int y = 0; y < indexFromProcNum(proc, I).y(); y++ )
        {
            si[proc].y() += Ni[procNumFromIndex(vector(0,y,0), I)].y();
        }
        for ( int z = 0; z < indexFromProcNum(proc, I).z(); z++ )
        {
            si[proc].z() += Ni[procNumFromIndex(vector(0,0,z), I)].z();
        }

        // Set starting indices of each processor, target decomp
        st[proc] = labelVector(0,0,0);
        for ( int x = 0; x < indexFromProcNum(proc, T).x(); x++ )
        {
            st[proc].x() += Nt[procNumFromIndex(vector(x,0,0), T)].x();
        }
        for ( int y = 0; y < indexFromProcNum(proc, T).y(); y++ )
        {
            st[proc].y() += Nt[procNumFromIndex(vector(0,y,0), T)].y();
        }
        for ( int z = 0; z < indexFromProcNum(proc, T).z(); z++ )
        {
            st[proc].z() += Nt[procNumFromIndex(vector(0,0,z), T)].z();
        }
    }

    List<labelVector> sendSize(Pstream::nProcs(),Zero);
    List<labelVector> sendStart(Pstream::nProcs(),Zero);

    List<labelVector> recvSize(Pstream::nProcs(),Zero);
    List<labelVector> recvStart(Pstream::nProcs(),Zero);

    // Loop over sending processors
    for (int i = 0; i < I.x(); i++)
    for (int j = 0; j < I.y(); j++)
    for (int k = 0; k < I.z(); k++)
    {
        labelVector ijk(i,j,k);
        label sendProcNum = procNumFromIndex(ijk,I);

        labelVector ei = si[sendProcNum] + Ni[sendProcNum];

        Info << nl;
        Info<< "Proc " << sendProcNum << " starting at " << si[sendProcNum]  << " and with size " << Ni[sendProcNum] << " sends to " << nl;

        // Loop over receiving processors
        for (int l = 0; l < T.x(); l++)
        for (int m = 0; m < T.y(); m++)
        for (int n = 0; n < T.z(); n++)
        {
            labelVector lmn(l,m,n);
            label recvProcNum = procNumFromIndex(lmn,T);
            // Start and end cell indices of receiving pocessor

            labelVector et = st[recvProcNum] + Nt[recvProcNum];

            // Check if send and recv processors overlap
            if
            (
                st[recvProcNum].x() < ei.x() &&
                st[recvProcNum].y() < ei.y() &&
                st[recvProcNum].z() < ei.z() &&
                et.x() > si[sendProcNum].x() &&
                et.y() > si[sendProcNum].y() &&
                et.z() > si[sendProcNum].z()
            )
            {
                // Overlap start and end cell indices

                labelVector s
                (
                    max(si[sendProcNum].x(), st[recvProcNum].x()),
                    max(si[sendProcNum].y(), st[recvProcNum].y()),
                    max(si[sendProcNum].z(), st[recvProcNum].z())
                );

                labelVector e
                (
                    min(ei.x(), et.x()),
                    min(ei.y(), et.y()),
                    min(ei.z(), et.z())
                );

                // Add send size and start to sendSize lists
                if (Pstream::myProcNo() == sendProcNum)
                {
                    sendSize[recvProcNum] = e - s;
                    sendStart[recvProcNum] = s;
                }

                // Add receive size and start to lists
                if (Pstream::myProcNo() == recvProcNum)
                {
                    recvSize[sendProcNum] = e - s;
                    recvStart[sendProcNum] = s;
                }

                Info<< "    Proc " << recvProcNum
                    << " -> Overlap starting at cell " << s
                    << " with size " << e - s << nl;
            }
        }
    }

    Info<< endl;

    // Prepare send buffer. Abuse the block class for this. Send a vector that
    // contains the global cell index. Sizes and displacements are in bytes
    // because we're sending as char.

    vectorBlock sendBuffer(Ni[Pstream::myProcNo()]);

    labelList sendCount(Pstream::nProcs());
    labelList sendDisplacement(Pstream::nProcs());

    // Input data, each block element contains its own global cell index
    vectorBlock inputData(Ni[Pstream::myProcNo()]);
    for ( int i = 0; i < Ni[Pstream::myProcNo()].x(); i++ )
    {
        for ( int j = 0; j < Ni[Pstream::myProcNo()].y(); j++ )
        {
            for ( int k = 0; k < Ni[Pstream::myProcNo()].z(); k++ )
            {
                inputData(i,j,k) = vector
                (
                    si[Pstream::myProcNo()].x() + i,
                    si[Pstream::myProcNo()].y() + j,
                    si[Pstream::myProcNo()].z() + k
                );
            }
        }
    }

    label cursor = 0;
    for (label proc = 0; proc < Pstream::nProcs(); proc++)
    {
        // Send displacement
        sendDisplacement[proc] = cursor*sizeof(vector);

        // Local coordinates of start point of the send buffer
        labelVector start
        (
            sendStart[proc].x() - si[Pstream::myProcNo()].x(),
            sendStart[proc].y() - si[Pstream::myProcNo()].y(),
            sendStart[proc].z() - si[Pstream::myProcNo()].z()
        );

        labelVector size = sendSize[proc];

        for (int i = start.x(); i < start.x() + size.x(); i++)
        for (int j = start.y(); j < start.y() + size.y(); j++)
        for (int k = start.z(); k < start.z() + size.z(); k++)
        {
            sendBuffer(cursor++) = inputData(i,j,k);
        }

        sendCount[proc] = cmptProduct(size)*sizeof(vector);
    }

    // Prepare receive buffer. Also abuse the block class for this.

    vectorBlock recvBuffer(Nt[Pstream::myProcNo()]);

    labelList recvCount(Pstream::nProcs());
    labelList recvDisplacement(Pstream::nProcs());

    cursor = 0;

    for (label proc = 0; proc < Pstream::nProcs(); proc++)
    {
        recvDisplacement[proc] = cursor*sizeof(vector);

        labelVector size = recvSize[proc];

        cursor += cmptProduct(size);

        recvCount[proc] = cmptProduct(size)*sizeof(vector);
    }

    // All-to-all using OpenFOAM's UPstream wrapper

    Info<< "Sending data with all-to-all ..." << endl;

    auto t1 = high_resolution_clock::now();

    for (int i = 0; i < repetitions; i++)
    {
        UPstream::allToAll
        (
            reinterpret_cast<char*>(sendBuffer.begin()),
            sendCount,
            sendDisplacement,
            reinterpret_cast<char*>(recvBuffer.begin()),
            recvCount,
            recvDisplacement,
            UPstream::worldComm
        );
    }

    auto t2 = high_resolution_clock::now();

    Info<< "Completed in "
        << duration_cast<milliseconds>(t2 - t1).count()
        << " ms" << endl;

    // Unpack data

    Info<< "Unpacking ..." << endl;

    vectorBlock data(Nt[Pstream::myProcNo()]);

    unpack(recvBuffer, recvStart, recvSize, data, st);

    Info<< "Data check ..." << endl;

    if (check(data, st[Pstream::myProcNo()]))
        Info<< "Data check passed" << endl;

    Info<< endl;


    // All-to-all using OpenFOAM's non-blocking functionality

    Info<< "Sending data with non-blocking transactions ..." << endl;

    data *= 0.0;

    t1 = high_resolution_clock::now();

    for (int i = 0; i < repetitions; i++)
    {
        for (int proc = 0; proc < Pstream::nProcs(); proc++)
        {
            if (recvCount[proc] > 0)
            {
                UIPstream::read
                (
                    Pstream::commsTypes::nonBlocking,
                    proc,
                    reinterpret_cast<char*>
                    (
                        &recvBuffer
                        (
                            recvDisplacement[proc]
                            / sizeof(vector)
                        )
                    ),
                    recvCount[proc],
                    0,
                    UPstream::worldComm
                );
            }

            if (sendCount[proc] > 0)
            {
                UOPstream::write
                (
                    Pstream::commsTypes::nonBlocking,
                    proc,
                    reinterpret_cast<char*>
                    (
                        &sendBuffer
                        (
                            sendDisplacement[proc]
                            / sizeof(vector)
                        )
                    ),
                    sendCount[proc],
                    0,
                    UPstream::worldComm
                );
            }
        }
    }

    UPstream::waitRequests();


    t2 = high_resolution_clock::now();

    Info<< "Completed in "
        << duration_cast<milliseconds>(t2 - t1).count()
        << " ms" << endl;

    // Unpack data

    Info<< "Unpacking ..." << endl;

    unpack(recvBuffer, recvStart, recvSize, data, st);

    Info<< "Data check ..." << endl;

    if (check(data, st[Pstream::myProcNo()]))
        Info<< "Data check passed" << endl;

    Info<< endl;


    // Send everything to processor 0
    Info << "Collecting data on processor 0..." << endl;

    vectorBlock totalRecvBuffer(N);

    labelList sc(Pstream::nProcs());
    labelList sd(Pstream::nProcs());
    labelList rc(Pstream::nProcs());
    labelList rd(Pstream::nProcs());

    for (label proc = 0; proc < Pstream::nProcs(); proc++)
    {
        sc[proc] = cmptProduct(Nt[Pstream::myProcNo()])*sizeof(vector);

        sd[proc] = 0;

        rc[proc] = cmptProduct(Nt[proc])*sizeof(vector);

        rd[proc] = 0;
        for (int p = 0; p < proc; p++)
        {
            rd[proc] += cmptProduct(Nt[p])*sizeof(vector);
        }
    }

    UPstream::allToAll
    (
        reinterpret_cast<char*>(data.begin()),
        sc,
        sd,
        reinterpret_cast<char*>(totalRecvBuffer.begin()),
        rc,
        rd,
        UPstream::worldComm
    );

    if ( ! Pstream::myProcNo() )
    {
        vectorBlock totalData(N);
        unpack(totalRecvBuffer, st, Nt, totalData, st);

        if(check(totalData, labelVector(0,0,0)))
            Info << "Total data check successful." << endl;
    }
}
