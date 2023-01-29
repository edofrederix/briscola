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
    vectorBlock& output
)
{
    label cursor = 0;

    labelVector N(output.shape());

    for (label proc = 0; proc < Pstream::nProcs(); proc++)
    {
        labelVector start
        (
            recvStart[proc].x() % N.x(),
            recvStart[proc].y() % N.y(),
            recvStart[proc].z() % N.z()
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
    if (N[dir] % I[dir] != 0)
    {
        FatalError
            << "Mesh not decomposable by initial decomposition vector"
            << endl;
        FatalError.exit();
    }

    for (int dir = 0; dir < 3; dir++)
    if (N[dir] % T[dir] != 0)
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

    const labelVector Ni = cmptDivide(N,I);
    const labelVector Nt = cmptDivide(N,T);

    List<labelVector> sendSize(Pstream::nProcs(),Zero);
    List<labelVector> sendStart(Pstream::nProcs(),Zero);

    List<labelVector> recvSize(Pstream::nProcs(),Zero);
    List<labelVector> recvStart(Pstream::nProcs(),Zero);

    for (int i = 0; i < I.x(); i++)
    for (int j = 0; j < I.y(); j++)
    for (int k = 0; k < I.z(); k++)
    {
        labelVector ijk(i,j,k);
        label sendProcNum = procNumFromIndex(ijk,I);

        // Start and end cell indices

        labelVector si(cmptMultiply(ijk,Ni));
        labelVector ei(cmptMultiply(ijk+unitXYZ,Ni));

        Info<< "Proc " << sendProcNum << " sends to " << nl;

        for (int l = si.x()/Nt.x(); l < (ei.x()-1)/Nt.x()+1; l++)
        for (int m = si.y()/Nt.y(); m < (ei.y()-1)/Nt.y()+1; m++)
        for (int n = si.z()/Nt.z(); n < (ei.z()-1)/Nt.z()+1; n++)
        {
            labelVector lmn(l,m,n);
            label recvProcNum = procNumFromIndex(lmn,T);

            // Start and end cell indices

            labelVector st(cmptMultiply(lmn,Nt));
            labelVector et(cmptMultiply(lmn+unitXYZ,Nt));

            // Overlap start and end cell indices

            labelVector s
            (
                max(si.x(), st.x()),
                max(si.y(), st.y()),
                max(si.z(), st.z())
            );

            labelVector e
            (
                min(ei.x(), et.x()),
                min(ei.y(), et.y()),
                min(ei.z(), et.z())
            );

            if (Pstream::myProcNo() == sendProcNum)
            {
                sendSize[recvProcNum] = e - s;
                sendStart[recvProcNum] = s;
            }

            if (Pstream::myProcNo() == recvProcNum)
            {
                recvSize[sendProcNum] = e - s;
                recvStart[sendProcNum] = s;
            }

            Info<< "    Proc " << recvProcNum
                << " starting at cell " << s
                << " with size " << e - s << nl;
        }
    }

    Info<< endl;

    // Prepare send buffer. Abuse the block class for this. Send a vector that
    // contains the global cell index. Sizes and displacements are in bytes
    // because we're sending as char.

    vectorBlock sendBuffer(Ni);

    labelList sendCount(Pstream::nProcs());
    labelList sendDisplacement(Pstream::nProcs());

    label cursor = 0;

    labelVector oi =
        cmptMultiply(indexFromProcNum(Pstream::myProcNo(), I), Ni);

    for (label proc = 0; proc < Pstream::nProcs(); proc++)
    {
        sendDisplacement[proc] = cursor*sizeof(vector);

        labelVector start
        (
            sendStart[proc].x() % Ni.x(),
            sendStart[proc].y() % Ni.y(),
            sendStart[proc].z() % Ni.z()
        );

        labelVector size = sendSize[proc];

        for (int i = start.x(); i < start.x() + size.x(); i++)
        for (int j = start.y(); j < start.y() + size.y(); j++)
        for (int k = start.z(); k < start.z() + size.z(); k++)
        {
            sendBuffer(cursor++) = vector
            (
                oi.x() + i,
                oi.y() + j,
                oi.z() + k
            );
        }

        sendCount[proc] = cmptProduct(size)*sizeof(vector);
    }

    // Prepare receive buffer. Also abuse the block class for this.

    vectorBlock recvBuffer(Nt);

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

    Info<< "Sending data wih all-to-all ..." << endl;

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

    vectorBlock data(Nt);

    unpack(recvBuffer, recvStart, recvSize, data);

    Info<< "Data check ..." << endl;

    labelVector ot =
        cmptMultiply(indexFromProcNum(Pstream::myProcNo(), T), Nt);

    if (check(data, ot))
        Info<< "Data check passed" << endl;

    Info<< endl;

    // All-to-all using OpenFOAM's non-blocking functionality

    Info<< "Sending data wih non-blocking transactions ..." << endl;

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

        UPstream::waitRequests();
    }

    t2 = high_resolution_clock::now();

    Info<< "Completed in "
        << duration_cast<milliseconds>(t2 - t1).count()
        << " ms" << endl;

    // Unpack data

    Info<< "Unpacking ..." << endl;

    unpack(recvBuffer, recvStart, recvSize, data);

    Info<< "Data check ..." << endl;

    if (check(data, ot))
        Info<< "Data check passed" << endl;

    Info<< endl;
}
