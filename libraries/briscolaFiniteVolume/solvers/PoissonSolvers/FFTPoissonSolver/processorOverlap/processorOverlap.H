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
    inline List<labelVector>& Ni() { return Ni_; }
    inline List<labelVector>& Si() { return Si_; }

    inline List<labelVector>& Nt() { return Nt_; }
    inline List<labelVector>& St() { return St_; }

    inline List<labelVector>& sendSize()  { return sendSize_;  }
    inline List<labelVector>& sendStart() { return sendStart_; }
    inline List<labelVector>& recvSize()  { return recvSize_;  }
    inline List<labelVector>& recvStart() { return recvStart_; }

    inline labelList& sendCount() { return sendCount_; }
    inline labelList& sendDisplacement() { return sendDisplacement_; }

    inline labelList& recvCount() { return recvCount_; }
    inline labelList& recvDisplacement() { return recvDisplacement_; }

};

}

}

}

#endif
