#include "decomposer.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

decomposer::decomposer
(
    const fvMesh& fvMsh,
    label decompType
)
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
    decompInit(decompType);
}

// Destructor

decomposer::~decomposer()
{}

// Initialize the initial and pencil decompositions

void decomposer::decompInit(label decompType)
{

    // Initial decomposition

    checkDecomp(I_);

    Ni_ = fvMsh_.msh().decomp().partSizePerProc();
    Si_ = fvMsh_.msh().decomp().globalStartPerProc();

    // Pencil decompositions

    switch (decompType)
    {
        case 0:
            X_ = Y_ = Z_ = I_;
            break;

        case 1:
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

            // Rotate X-pencils to Y-pencils
            Y_ = vector(X_.y(), 1, X_.z());

            // Rotate Y-pencils to Z-pencils
            Z_ = vector(Y_.x(), Y_.z(), 1);

            break;

        case 2:
            X_ = I_;

            // Rotate X-pencils to Y-pencils
            Y_ = vector(X_.y(), 1, X_.z());

            // Rotate Y-pencils to Z-pencils
            Z_ = vector(Y_.x(), Y_.z(), 1);

            break;

        case 3:
            Y_ = I_;

            // Rotate Y-pencils to Z-pencils
            Z_ = vector(Y_.x(), Y_.z(), 1);

            // Rotate Z-pencils to X-pencils
            X_ = vector(1, Z_.x(), Z_.y());

            break;

        case 4:
            Z_ = I_;

            // Rotate Z-pencils to X-pencils
            X_ = vector(1, Z_.x(), Z_.y());

            // Rotate X-pencils to Y-pencils
            Y_ = vector(X_.y(), 1, X_.z());

            break;

        case 5:
            X_ = Y_ = I_;

            // Rotate XY-slab to XZ-slab
            Z_ = vector(1, Y_.z(), 1);

            break;

        case 6:
            X_ = Z_ = I_;

            // Rotate XZ-slab to XY-slab
            Y_ = vector(1, 1, X_.y());

            break;

        case 7:
            Y_ = Z_ = I_;

            // Rotate YZ-slab to XZ-slab
            X_ = vector(1, Z_.x(), 1);

            break;
    }

    // Check that all decompositions are valid
    checkDecomp(X_);
    checkDecomp(Y_);
    checkDecomp(Z_);

    // Set dimensions and origins

    Nx_ = procDims(X_);
    Sx_ = procOrig(X_);

    Ny_ = procDims(Y_);
    Sy_ = procOrig(Y_);

    Nz_ = procDims(Z_);
    Sz_ = procOrig(Z_);
}

void decomposer::checkDecomp(labelVector D)
{
    if (cmptProduct(D) != Pstream::nProcs())
    {
        FatalError
            << "Decomposition " << D << " does not match number of processors."
            << endl;
        FatalError.exit();
    }

    for (int dir = 0; dir < 3; dir++)
    if (N_[dir] < D[dir])
    {
        FatalError
            << "Mesh not decomposable by decomposition vector " << D
            << endl;
        FatalError.exit();
    }
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
// Only works when I != T

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

            // If remainder > 0, the data distribution is uneven

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
    labelVector T,
    string recvMajorOrder
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
            if (recvMajorOrder == "x")
            {
                output(k*N.x()*N.y() + j*N.x() + i) = buffer(cursor++);
            }
            else if (recvMajorOrder == "y")
            {
                output(i*N.y()*N.z() + k*N.y() + j) = buffer(cursor++);
            }
            else
            {
                output(i,j,k) = buffer(cursor++);
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
    labelVector T,
    string sendMajorOrder,
    string recvMajorOrder
)
{
    if (I == T)
    {
        if (sendMajorOrder == recvMajorOrder)
        {
            dst = src;
        }
        // Only change the major order of the array
        else
        {
            labelVector Nd = src.shape();

            scalarBlock buffer(Nd, Zero);

            if (sendMajorOrder == "x")
            {
                for (int i = 0; i < Nd.x(); i++)
                for (int j = 0; j < Nd.y(); j++)
                for (int k = 0; k < Nd.z(); k++)
                {
                    buffer(i,j,k) = src
                    (
                        k * Nd.x()
                            * Nd.y()
                        + j * Nd.x()
                        + i
                    );
                }
            }
            else if (sendMajorOrder == "y")
            {
                for (int i = 0; i < Nd.x(); i++)
                for (int j = 0; j < Nd.y(); j++)
                for (int k = 0; k < Nd.z(); k++)
                {
                    buffer(i,j,k) = src
                    (
                        i * Nd.y()
                            * Nd.z()
                        + k * Nd.y()
                        + j
                    );
                }
            }
            else
            {
                buffer = src;
            }

            // Most arguments of unpack can be set
            // to zero due to complete overlap

            List<labelVector> recvStart(Pstream::nProcs(), Zero);
            List<labelVector> recvSize(Pstream::nProcs(), Zero);
            List<labelVector> St(Pstream::nProcs(), Zero);

            // Only the recvSize of the current processor is needed

            recvSize[Pstream::myProcNo()] = src.shape();

            unpack(buffer, recvStart, recvSize, dst, St, T, recvMajorOrder);
        }

    }
    else
    {
        // Make decomposition pair key

        word IT = makeKey(I,T);

        // Add processorOverlap to table if it isn't there already

        if (!processorOverlap_.found(IT))
        {
            processorOverlap overlap(*this, I, T);
            processorOverlap_.insert(IT, overlap);
        }

        // Processor overlap variables

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
                // x-pencil is stored as k>j>i so that elements
                // in a line in x are close to each other in memory
                if (sendMajorOrder == "x")
                {
                    sendBuffer(cursor++) = src
                    (
                          k * Ni[Pstream::myProcNo()].x()
                            * Ni[Pstream::myProcNo()].y()
                        + j * Ni[Pstream::myProcNo()].x()
                        + i
                    );
                }
                // y-pencil is stored as i>k>j
                else if (sendMajorOrder == "y")
                {
                    sendBuffer(cursor++) = src
                    (
                          i * Ni[Pstream::myProcNo()].y()
                            * Ni[Pstream::myProcNo()].z()
                        + k * Ni[Pstream::myProcNo()].y()
                        + j
                    );
                }
                // Other decompositions are stored as i>j>k
                else
                {
                    sendBuffer(cursor++) = src(i,j,k);
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

        unpack(recvBuffer, recvStart, recvSize, dst, St, T, recvMajorOrder);
    }
}

} // end namespace fv

} // end namespace briscola

} // end namespace Foam