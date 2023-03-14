#ifndef vertexVector_H
#define vertexVector_H

#include "vector.H"
#include "VertexSpace.H"
#include "contiguous.H"

namespace Foam
{

typedef VertexSpace<vector> vertexVector;

template<>
inline bool contiguous<vertexVector>() { return true; }

}

#endif