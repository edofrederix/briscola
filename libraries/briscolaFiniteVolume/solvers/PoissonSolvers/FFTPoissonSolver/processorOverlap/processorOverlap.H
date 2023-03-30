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

class decomposer;

class processorOverlap
{

    decomposer& decomp_;

    labelVector I_;
    labelVector T_;

    List<labelVector> Ni_;
    List<labelVector> Si_;

    List<labelVector> Nt_;
    List<labelVector> St_;

    List<labelVector> sendSize_;
    List<labelVector> sendStart_;

    List<labelVector> recvSize_;
    List<labelVector> recvStart_;

    labelList sendCount_;
    labelList sendDisplacement_;

    labelList recvCount_;
    labelList recvDisplacement_;

    // Processor overlap
    void computeOverlap();

public:

    // Construct from dictionary
    processorOverlap
    (
        decomposer& d,
        labelVector I,
        labelVector T
    );

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

} // end namespace fv

} // end namespace briscola

} // end namespace Foam

#endif
