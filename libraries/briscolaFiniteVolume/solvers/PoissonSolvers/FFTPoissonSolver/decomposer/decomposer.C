#include "decomposer.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

decomposer::decomposer(const fvMesh& fvMsh)
:
    fvMsh_(fvMsh),
    N_(fvMsh_.msh().cast<rectilinearMesh>().N()),
    I_(fvMsh_.msh().decomp().map().legend()[Pstream::nProcs() - 1] + unitXYZ)
{
    Ni_ = List<labelVector>(Pstream::nProcs(), Zero);
    Nx_ = List<labelVector>(Pstream::nProcs(), Zero);
    Ny_ = List<labelVector>(Pstream::nProcs(), Zero);
    Nz_ = List<labelVector>(Pstream::nProcs(), Zero);

    si_ = List<labelVector>(Pstream::nProcs(), Zero);
    sx_ = List<labelVector>(Pstream::nProcs(), Zero);
    sy_ = List<labelVector>(Pstream::nProcs(), Zero);
    sz_ = List<labelVector>(Pstream::nProcs(), Zero);

    decompInit();
}

// Destructor

decomposer::~decomposer()
{}

// Initialize the initial and pencil decompositions

void decomposer::decompInit()
{
    // Some checks

    if (cmptProduct(I_) != Pstream::nProcs())
    {
        FatalError
            << "Initial decomposition does not match the number of processors"
            << endl;
        FatalError.exit();
    }

    for (int dir = 0; dir < 3; dir++)
    if (N_[dir] < I_[dir])
    {
        FatalError
            << "Mesh not decomposable by initial decomposition vector"
            << endl;
        FatalError.exit();
    }

    // Print initial processor topology

    Info<< "Initial topology: " << nl;

    for (label proc = 0; proc < Pstream::nProcs(); proc++)
        Info<< "    Proc " << proc << " at "
            << fvMsh_.msh().decomp().map().legend()[proc] << endl;

    Info<< endl;

    // Initial decomposition

    Ni_ = fvMsh_.msh().decomp().partSizePerProc();
    si_ = fvMsh_.msh().decomp().globalStartPerProc();

    // Pencil decompositions: check if existing decompositions may be reused

    if (I_.x() > 1)
    {
        decompX();
    }
    else
    {
        X_ = I_;
        Nx_ = Ni_;
        sx_ = si_;
    }

    if (X_.y() > 1)
    {
        decompY();
    }
    else
    {
        Y_ = X_;
        Ny_ = Nx_;
        sy_ = sx_;
    }

    if (Y_.z() > 1)
    {
        decompZ();
    }
    else
    {
        Z_ = Y_;
        Nz_ = Ny_;
        sz_ = sy_;
    }
}

// Initialize X-pencil decomposition

void decomposer::decompX()
{
    // Factorize number of processors for x-pencil deocmposition

    X_ = vector
    (
        1,
        ceil( Foam::sqrt( static_cast<scalar>( Pstream::nProcs() ) ) ),
        0
    );

    while (Pstream::nProcs() % X_.y() != 0)
    {
        X_.y() -= 1;
    }

    X_.z() = Pstream::nProcs() / X_.y();

    // Check decomposition

    for (int dir = 0; dir < 3; dir++)
    if (N_[dir] < X_[dir])
    {
        FatalError
            << "X-pencil decomposition failed."
            << endl;
        FatalError.exit();
    }

    // Set dimensions and origins

    Nx_ = procDims(X_);
    sx_ = procOrig(X_);

    Info<< "X-pencil topology: " << nl;

    for (label proc = 0; proc < Pstream::nProcs(); proc++)
        Info<< "    Proc " << proc << " at "
            << indexFromProcNum(proc, X_) << endl;

    Info<< endl;

}

// Initialize Y-pencil decomposition

void decomposer::decompY()
{
    // Rotate X-pencils to Y-pencils

    Y_ = vector
    (
        X_.y(),
        1,
        X_.z()
    );

    // Check decomposition

    for (int dir = 0; dir < 3; dir++)
    if (N_[dir] < Y_[dir])
    {
        FatalError
            << "Y-pencil decomposition failed."
            << endl;
        FatalError.exit();
    }

    // Set dimensions and origins

    Ny_ = procDims(Y_);
    sy_ = procOrig(Y_);

    Info<< "Y-pencil topology: " << nl;

    for (label proc = 0; proc < Pstream::nProcs(); proc++)
        Info<< "    Proc " << proc << " at "
            << indexFromProcNum(proc, Y_) << endl;

    Info<< endl;

}

// Initialize Z-pencil decomposition

void decomposer::decompZ()
{
    // Rotate Y-pencils to Z-pencils

    Z_ = vector
    (
        Y_.x(),
        Y_.z(),
        1
    );

    // Check decomposition

    for (int dir = 0; dir < 3; dir++)
    if (N_[dir] < Z_[dir])
    {
        FatalError
            << "Z-pencil decomposition failed."
            << endl;
        FatalError.exit();
    }

    // Set dimensions and origins

    Nz_ = procDims(Z_);
    sz_ = procOrig(Z_);

    Info<< "Z-pencil topology: " << nl;

    for (label proc = 0; proc < Pstream::nProcs(); proc++)
        Info<< "    Proc " << proc << " at "
            << indexFromProcNum(proc, Z_) << endl;

    Info<< endl;

}

// Return the processor number given an index and decomposition

label decomposer::procNumFromIndex
(
    const labelVector ijk,
    const labelVector D
)
{
    if (D == I_)
    {
        return fvMsh_.msh().decomp().map()(ijk);
    }
    else
    {
        return ijk.x()*D.y()*D.z() + ijk.y()*D.z() + ijk.z();
    }
}

// Return the index given a processor number and decomposition

labelVector decomposer::indexFromProcNum
(
    const label num,
    const labelVector D
)
{
    if (D == I_)
    {
        return fvMsh_.msh().decomp().map().legend()[num];
    }
    else
    {
        label i = num/(D.y()*D.z());
        label j = (num-i*D.y()*D.z())/D.z();
        label k = (num-i*D.y()*D.z()-j*D.z());

        return labelVector(i,j,k);
    }
}

// Return the list of processor dimensions for a given decomposition

List<labelVector> decomposer::procDims(labelVector D)
{
    // Some checks

    if (cmptProduct(D) != Pstream::nProcs())
    {
        FatalError
            << "Initial decomposition does not match the number of processors"
            << endl;
        FatalError.exit();
    }

    for (int dir = 0; dir < 3; dir++)
    if (N_[dir] < D[dir])
    {
        FatalError
            << "Mesh not decomposable by initial decomposition vector"
            << endl;
        FatalError.exit();
    }

    if (D == I_)
    {
        return Ni_;
    }
    else
    {
        List<labelVector> Nd(Pstream::nProcs(), Zero);

        for ( int proc = 0; proc < Pstream::nProcs(); proc++ )
        {
            Nd[proc] = cmptDivide(N_,D);

            // If the data distribution is uneven, remainder > 0

            labelVector remainder = N_ - cmptMultiply(Nd[proc], D);

            Nd[proc] += labelVector
                (
                    indexFromProcNum(proc,D).x() < remainder.x(),
                    indexFromProcNum(proc,D).y() < remainder.y(),
                    indexFromProcNum(proc,D).z() < remainder.z()
                );
        }

        return Nd;
    }
}

// Return list of processor origin indices for a given decomposition

List<labelVector> decomposer::procOrig(labelVector D)
{
    if (D == I_)
    {
        return si_;
    }
    else
    {
        List<labelVector> Nd(procDims(D));
        List<labelVector> sd = List<labelVector>(Pstream::nProcs(), Zero);

        for ( int proc = 0; proc < Pstream::nProcs(); proc++ )
        {

            // Set starting indices of each processor

            for ( int x = 0; x < indexFromProcNum(proc, D).x(); x++ )
                sd[proc].x() +=
                    Nd[procNumFromIndex(labelVector(x,0,0), D)].x();

            for ( int y = 0; y < indexFromProcNum(proc, D).y(); y++ )
                sd[proc].y() +=
                    Nd[procNumFromIndex(labelVector(0,y,0), D)].y();

            for ( int z = 0; z < indexFromProcNum(proc, D).z(); z++ )
                sd[proc].z() +=
                    Nd[procNumFromIndex(labelVector(0,0,z), D)].z();

        }

        return sd;
    }
}

// Unpack a receive buffer

void decomposer::unpack
(
    const scalarBlock& buffer,
    const List<labelVector>& recvStart,
    const List<labelVector>& recvSize,
    scalarBlock& output,
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

// Transpose data from a source to a destination block,
// from an initial decomposition to a target decomposition

void decomposer::transpose
(
    const scalarBlock& src,
    scalarBlock& dst,
    labelVector I,
    labelVector T
)
{
    if (I == T)
    {
        Info << "Initial and target decompositions identical." << nl;

        dst = src;
    }
    else
    {
        List<labelVector> Ni = procDims(I);
        List<labelVector> si = procOrig(I);

        List<labelVector> Nt = procDims(T);
        List<labelVector> st = procOrig(T);

        // Compute processor overlap from initial to target,
        // defined by the global start cell index and block size

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
                    st[recvProcNum].x() < ei.x()
                && st[recvProcNum].y() < ei.y()
                && st[recvProcNum].z() < ei.z()
                && et.x() > si[sendProcNum].x()
                && et.y() > si[sendProcNum].y()
                && et.z() > si[sendProcNum].z()
                )
                {
                    // Overlap start and end cell indices

                    labelVector s
                    (
                        std::max(si[sendProcNum].x(), st[recvProcNum].x()),
                        std::max(si[sendProcNum].y(), st[recvProcNum].y()),
                        std::max(si[sendProcNum].z(), st[recvProcNum].z())
                    );

                    labelVector e
                    (
                        std::min(ei.x(), et.x()),
                        std::min(ei.y(), et.y()),
                        std::min(ei.z(), et.z())
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
                }
            }

            // Info << endl;
        }

        // Prepare send buffer. Abuse the block class for this. Send a vector
        // that contains the global cell index. Sizes and displacements are
        // in bytes because we're sending as char.

        scalarBlock sendBuffer(Ni[Pstream::myProcNo()]);

        labelList sendCount(Pstream::nProcs());
        labelList sendDisplacement(Pstream::nProcs());

        label cursor = 0;

        for (label proc = 0; proc < Pstream::nProcs(); proc++)
        {
            // Send displacement

            sendDisplacement[proc] = cursor*sizeof(scalar);

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
                sendBuffer(cursor++) = src(i,j,k);
            }

            sendCount[proc] = cmptProduct(size)*sizeof(scalar);
        }

        // Prepare receive buffer. Also abuse the block class for this.

        scalarBlock recvBuffer(Nt[Pstream::myProcNo()]);

        labelList recvCount(Pstream::nProcs());
        labelList recvDisplacement(Pstream::nProcs());

        cursor = 0;

        for (label proc = 0; proc < Pstream::nProcs(); proc++)
        {
            recvDisplacement[proc] = cursor*sizeof(scalar);

            labelVector size = recvSize[proc];

            cursor += cmptProduct(size);

            recvCount[proc] = cmptProduct(size)*sizeof(scalar);
        }

        // All-to-all using OpenFOAM's UPstream wrapper

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

        unpack(recvBuffer, recvStart, recvSize, dst, st);
    }
}

// Transform from i>j>k to i>k>j

void decomposer::yTransFwd(scalarBlock& yData)
{
    scalarBlock yCopy = yData;

    labelVector N(yData.shape());

    label cursor = 0;
    for (int i = 0; i < N.x(); i++)
    {
        for (int k = 0; k < N.z(); k++)
        {
            for (int j = 0; j < N.y(); j++)
            {
                yData(cursor++) = yCopy(i,j,k);
            }
        }
    }
}

// Transform from i>k>j to i>j>k

void decomposer::yTransBwd(scalarBlock& yData)
{
    scalarBlock yCopy = yData;

    labelVector N(yData.shape());

    label cursor = 0;
    for (int i = 0; i < N.x(); i++)
    {
        for (int k = 0; k < N.z(); k++)
        {
            for (int j = 0; j < N.y(); j++)
            {
                yData(i,j,k) = yCopy(cursor++);
            }
        }
    }

}

} // end namespace fv

} // end namespace briscola

} // end namespace Foam