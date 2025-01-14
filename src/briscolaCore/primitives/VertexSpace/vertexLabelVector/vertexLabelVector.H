#ifndef vertexLabelVector_H
#define vertexLabelVector_H

#include "labelVector.H"
#include "VertexSpace.H"
#include "contiguous.H"

namespace Foam
{

typedef VertexSpace<labelVector> vertexLabelVector;

template<>
inline bool contiguous<vertexLabelVector>() { return true; }

}

#endif