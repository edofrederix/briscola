#include "processorOverlap.H"
#include "pencilDecomposer.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace FFT
{

void processorOverlap::computeOverlap()
{
    Ni_ = decomp_.procDims(I_);
    Si_ = decomp_.procOrig(I_);

    Nt_ = decomp_.procDims(T_);
    St_ = decomp_.procOrig(T_);

    // Loop over sending processors

    for (int i = 0; i < I_.x(); i++)
    for (int j = 0; j < I_.y(); j++)
    for (int k = 0; k < I_.z(); k++)
    {
        labelVector ijk(i,j,k);
        label sendProcNum = decomp_.procNumFromIndex(ijk,I_);

        labelVector ei = Si_[sendProcNum] + Ni_[sendProcNum];

        // Loop over receiving processors

        for (int l = 0; l < T_.x(); l++)
        for (int m = 0; m < T_.y(); m++)
        for (int n = 0; n < T_.z(); n++)
        {
            labelVector lmn(l,m,n);
            label recvProcNum = decomp_.procNumFromIndex(lmn,T_);

            // Start and end cell indices of receiving pocessor

            labelVector et = St_[recvProcNum] + Nt_[recvProcNum];

            // Check if send and recv processors overlap

            if
            (
                   St_[recvProcNum].x() < ei.x()
                && St_[recvProcNum].y() < ei.y()
                && St_[recvProcNum].z() < ei.z()
                && et.x() > Si_[sendProcNum].x()
                && et.y() > Si_[sendProcNum].y()
                && et.z() > Si_[sendProcNum].z()
            )
            {
                // Overlap start and end cell indices

                labelVector s
                (
                    std::max(Si_[sendProcNum].x(), St_[recvProcNum].x()),
                    std::max(Si_[sendProcNum].y(), St_[recvProcNum].y()),
                    std::max(Si_[sendProcNum].z(), St_[recvProcNum].z())
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
                    sendSize_[recvProcNum] = e - s;
                    sendStart_[recvProcNum] = s;
                }

                // Add receive size and start to lists

                if (Pstream::myProcNo() == recvProcNum)
                {
                    recvSize_[sendProcNum] = e - s;
                    recvStart_[sendProcNum] = s;
                }
            }
        }
    }

    // Sizes and displacements are in
    // bytes because we're sending as char

    label cursor = 0;

    for (label proc = 0; proc < Pstream::nProcs(); proc++)
    {
        sendDisplacement_[proc] = cursor*sizeof(scalar);

        labelVector size = sendSize_[proc];

        sendCount_[proc] = cmptProduct(size)*sizeof(scalar);

        cursor += cmptProduct(size);
    }

    cursor = 0;

    for (label proc = 0; proc < Pstream::nProcs(); proc++)
    {
        recvDisplacement_[proc] = cursor*sizeof(scalar);

        labelVector size = recvSize_[proc];

        recvCount_[proc] = cmptProduct(size)*sizeof(scalar);

        cursor += cmptProduct(size);
    }
}

// Constructor

processorOverlap::processorOverlap
(
    pencilDecomposer& d,
    labelVector I,
    labelVector T
)
:
    decomp_(d),
    I_(I),
    T_(T),
    sendSize_(Pstream::nProcs(),Zero),
    sendStart_(Pstream::nProcs(),Zero),
    recvSize_(Pstream::nProcs(),Zero),
    recvStart_(Pstream::nProcs(),Zero),
    sendCount_(Pstream::nProcs()),
    sendDisplacement_(Pstream::nProcs()),
    recvCount_(Pstream::nProcs()),
    recvDisplacement_(Pstream::nProcs())
{
    computeOverlap();
}

// Copy constructor

processorOverlap::processorOverlap
(
    const processorOverlap& po
)
:
    decomp_(po.decomp_),
    I_(po.I_),
    T_(po.T_),
    Ni_(po.Ni_),
    Si_(po.Si_),
    Nt_(po.Nt_),
    St_(po.St_),
    sendSize_(po.sendSize_),
    sendStart_(po.sendStart_),
    recvSize_(po.recvSize_),
    recvStart_(po.recvStart_),
    sendCount_(po.sendCount_),
    sendDisplacement_(po.sendDisplacement_),
    recvCount_(po.recvCount_),
    recvDisplacement_(po.recvDisplacement_)
{}

// Destructor

processorOverlap::~processorOverlap()
{}

}

}

}

}
