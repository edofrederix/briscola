#ifndef lowerFaceVector_H
#define lowerFaceVector_H

#include "vector.H"
#include "LowerFaceSpace.H"
#include "contiguous.H"

namespace Foam
{

typedef LowerFaceSpace<vector> lowerFaceVector;

template<>
inline bool contiguous<lowerFaceVector>() { return true; }

}

#endif