#ifndef edgeVector_H
#define edgeVector_H

#include "vector.H"
#include "EdgeSpace.H"
#include "contiguous.H"

namespace Foam
{

typedef EdgeSpace<vector> edgeVector;

template<>
inline bool contiguous<edgeVector>() { return true; }

}

#endif