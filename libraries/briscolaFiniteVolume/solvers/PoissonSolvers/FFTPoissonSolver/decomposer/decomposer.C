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
    I_(fvMsh_.msh().decomp().map().legend()[Pstream::nProcs() - 1] + unitXYZ),
    Ni_(Pstream::nProcs(), Zero),
    Nx_(Pstream::nProcs(), Zero),
    Ny_(Pstream::nProcs(), Zero),
    Nz_(Pstream::nProcs(), Zero),
    Si_(Pstream::nProcs(), Zero),
    Sx_(Pstream::nProcs(), Zero),
    Sy_(Pstream::nProcs(), Zero),
    Sz_(Pstream::nProcs(), Zero)
{
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

    // Initial decomposition

    Ni_ = fvMsh_.msh().decomp().partSizePerProc();
    Si_ = fvMsh_.msh().decomp().globalStartPerProc();

    // Pencil decompositions: check if existing decompositions may be reused

    if (I_.x() > 1)
    {
        decompX();
    }
    else
    {
        X_ = I_;
        Nx_ = Ni_;
        Sx_ = Si_;
    }

    if (X_.y() > 1)
    {
        decompY();
    }
    else
    {
        Y_ = X_;
        Ny_ = Nx_;
        Sy_ = Sx_;
    }

    if (Y_.z() > 1)
    {
        decompZ();
    }
    else
    {
        Z_ = Y_;
        Nz_ = Ny_;
        Sz_ = Sy_;
    }
}

// Initialize X-pencil decomposition

void decomposer::decompX()
{
    // Factorize number of processors for x-pencil deocmposition

    X_ = vector
    (
        1,
        ceil(Foam::sqrt(scalar(Pstream::nProcs()))),
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
    Sx_ = procOrig(X_);
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
    Sy_ = procOrig(Y_);
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
    Sz_ = procOrig(Z_);
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

// Make key

word decomposer::makeKey(labelVector I, labelVector T)
{
    word key = "";

    if (I == I_)
    {
        key += 'I';
    }
    else if (I == X_)
    {
        key += 'X';
    }
    else if (I == Y_)
    {
        key += 'Y';
    }
    else if (I == Z_)
    {
        key += 'Z';
    }

    if (T ==  I_)
    {
        key += 'I';
    }
    else if (T ==  X_)
    {
        key += 'X';
    }
    else if (T ==  Y_)
    {
        key += 'Y';
    }
    else if (T ==  Z_)
    {
        key += 'Z';
    }

    return key;
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

        for (int proc = 0; proc < Pstream::nProcs(); proc++)
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
        return Si_;
    }
    else
    {
        List<labelVector> Nd(procDims(D));
        List<labelVector> Sd = List<labelVector>(Pstream::nProcs(), Zero);

        for (int proc = 0; proc < Pstream::nProcs(); proc++)
        {

            // Set starting indices of each processor

            for (int x = 0; x < indexFromProcNum(proc, D).x(); x++)
                Sd[proc].x() +=
                    Nd[procNumFromIndex(labelVector(x,0,0), D)].x();

            for (int y = 0; y < indexFromProcNum(proc, D).y(); y++)
                Sd[proc].y() +=
                    Nd[procNumFromIndex(labelVector(0,y,0), D)].y();

            for (int z = 0; z < indexFromProcNum(proc, D).z(); z++)
                Sd[proc].z() +=
                    Nd[procNumFromIndex(labelVector(0,0,z), D)].z();
        }

        return Sd;
    }
}

// Unpack a receive buffer

void decomposer::unpack
(
    const scalarBlock& buffer,
    const List<labelVector>& recvStart,
    const List<labelVector>& recvSize,
    scalarBlock& output,
    List<labelVector>& startIndex,
    bool yPencil
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
            if (!yPencil)
            {
                output(i,j,k) = buffer(cursor++);
            }
            // y-pencil decomposition is stored as i>k>j
            // due to the data format requirements of FFTW
            else
            {
                output(i*N.y()*N.z() + k*N.y() + j) = buffer(cursor++);
            }
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
        dst = src;
    }
    else
    {
        word IT = makeKey(I,T);

        if (!processorOverlap_.found(IT))
        {
            processorOverlap overlap(*this, I, T);
            processorOverlap_.insert(IT, overlap);
        }

        List<labelVector>& Ni = processorOverlap_[IT].Ni();
        List<labelVector>& Si = processorOverlap_[IT].Si();

        List<labelVector>& Nt = processorOverlap_[IT].Nt();
        List<labelVector>& St = processorOverlap_[IT].St();

        List<labelVector>& sendSize = processorOverlap_[IT].sendSize();
        List<labelVector>& sendStart = processorOverlap_[IT].sendStart();
        List<labelVector>& recvSize = processorOverlap_[IT].recvSize();
        List<labelVector>& recvStart = processorOverlap_[IT].recvStart();

        labelList& sendCount = processorOverlap_[IT].sendCount();
        labelList& sendDisplacement = processorOverlap_[IT].sendDisplacement();

        labelList& recvCount = processorOverlap_[IT].recvCount();
        labelList& recvDisplacement = processorOverlap_[IT].recvDisplacement();

        // Prepare send buffer

        scalarBlock sendBuffer(Ni[Pstream::myProcNo()]);

        label cursor = 0;

        for (label proc = 0; proc < Pstream::nProcs(); proc++)
        {
            // Local coordinates of start point of the send buffer

            labelVector start
            (
                sendStart[proc].x() - Si[Pstream::myProcNo()].x(),
                sendStart[proc].y() - Si[Pstream::myProcNo()].y(),
                sendStart[proc].z() - Si[Pstream::myProcNo()].z()
            );

            labelVector size = sendSize[proc];

            for (int i = start.x(); i < start.x() + size.x(); i++)
            for (int j = start.y(); j < start.y() + size.y(); j++)
            for (int k = start.z(); k < start.z() + size.z(); k++)
            {
                if (I != Y_)
                {
                    sendBuffer(cursor++) = src(i,j,k);
                }
                // y-pencil decomposition is stored as i>k>j
                // due to the data format requirements of FFTW
                else
                {
                    sendBuffer(cursor++) = src
                    (
                          i * Ni[Pstream::myProcNo()].y()
                            * Ni[Pstream::myProcNo()].z()
                        + k * Ni[Pstream::myProcNo()].y()
                        + j
                    );
                }
            }
        }

        // Prepare receive buffer

        scalarBlock recvBuffer(Nt[Pstream::myProcNo()]);

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

        if (T != Y_)
        {
            unpack(recvBuffer, recvStart, recvSize, dst, St);
        }
        // y-pencil decomposition is stored as i>k>j
        // due to the data format requirements of FFTW
        else
        {
            unpack(recvBuffer, recvStart, recvSize, dst, St, true);
        }
    }
}

} // end namespace fv

} // end namespace briscola

} // end namespace Foam