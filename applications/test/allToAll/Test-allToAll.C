#include "arguments.H"

#include "vector.H"
#include "block.H"

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

int main(int argc, char *argv[])
{
    arguments::addBoolOption("parallel", "run in parallel");

    arguments::validArgs.append("mesh size");
    arguments::validArgs.append("initial decomposition");
    arguments::validArgs.append("target decomposition");

    arguments args(argc, argv);

    labelVector N(args.argRead<labelVector>(1));
    labelVector I(args.argRead<labelVector>(2));
    labelVector T(args.argRead<labelVector>(3));

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

    Info<< "Sending data" << endl << endl;

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

    // Unpack data

    vectorBlock data(Nt);

    cursor = 0;

    for (label proc = 0; proc < Pstream::nProcs(); proc++)
    {
        labelVector start
        (
            recvStart[proc].x() % Nt.x(),
            recvStart[proc].y() % Nt.y(),
            recvStart[proc].z() % Nt.z()
        );

        labelVector size = recvSize[proc];

        for (int i = start.x(); i < start.x() + size.x(); i++)
        for (int j = start.y(); j < start.y() + size.y(); j++)
        for (int k = start.z(); k < start.z() + size.z(); k++)
        {
            data(i,j,k) = recvBuffer(cursor++);
        }
    }

    // Test by checking if the received vectors match the cell indices

    labelVector ot =
        cmptMultiply(indexFromProcNum(Pstream::myProcNo(), T), Nt);

    forAllBlock(data, i, j, k)
    {
        labelVector test = ot + labelVector(i,j,k);

        if (labelVector(data(i,j,k)) != test)
            FatalError
                << "Test failed. Expected " << test
                << " but received " << data(i,j,k) << endl
                << "At index " << labelVector(i,j,k) << endl
                << abort(FatalError);
    }

    Info<< "Data test succeeded" << endl << endl;
}
