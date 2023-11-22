#ifndef edgeLabelVector_H
#define edgeLabelVector_H

#include "labelVector.H"
#include "EdgeSpace.H"
#include "contiguous.H"

namespace Foam
{

typedef EdgeSpace<labelVector> edgeLabelVector;

template<>
inline bool contiguous<edgeLabelVector>() { return true; }

}

#endif
