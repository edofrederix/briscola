#ifndef processorOverlap_H
#define processorOverlap_H

#include "List.H"
#include "scalarList.H"
#include "vector.H"
#include "block.H"
#include "fvMesh.H"
#include "rectilinearMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace FFT
{

class pencilDecomposer;

class processorOverlap
{
    // Decomposer object
    pencilDecomposer& decomp_;

    // Initial and target decomposition vectors
    labelVector I_;
    labelVector T_;

    // List of processor part dimensions and starting indices
    // for initial decomposition and for target decomposition
    List<labelVector> Ni_;
    List<labelVector> Si_;

    List<labelVector> Nt_;
    List<labelVector> St_;

    // Size and start index of block part to be sent to each processor
    List<labelVector> sendSize_;
    List<labelVector> sendStart_;

    // Size and start index of block part to be received from each processor
    List<labelVector> recvSize_;
    List<labelVector> recvStart_;

    // Send size and displacement in bytes
    labelList sendCount_;
    labelList sendDisplacement_;

    // Receive count and displacement in bytes
    labelList recvCount_;
    labelList recvDisplacement_;

    // Compute processor overlaps
    void computeOverlap();

public:

    // Construct from decomposition parir
    processorOverlap
    (
        pencilDecomposer& d,
        labelVector I,
        labelVector T
    );

    // Copy constructor
    processorOverlap(const processorOverlap&);

    // Destructor
    virtual ~processorOverlap();

    // Inline return fucntions
    inline const List<labelVector>& Ni() const { return Ni_; }
    inline const List<labelVector>& Si() const { return Si_; }

    inline const List<labelVector>& Nt() const { return Nt_; }
    inline const List<labelVector>& St() const { return St_; }

    inline const List<labelVector>& sendSize() const { return sendSize_;  }
    inline const List<labelVector>& sendStart() const { return sendStart_; }
    inline const List<labelVector>& recvSize() const { return recvSize_;  }
    inline const List<labelVector>& recvStart() const { return recvStart_; }

    inline const labelList& sendCount() const { return sendCount_; }
    inline const labelList& sendDisplacement() const
    {
        return sendDisplacement_;
    }

    inline const labelList& recvCount() const { return recvCount_; }
    inline const labelList& recvDisplacement() const
    {
        return recvDisplacement_;
    }

};

}

}

}

}

#endif
