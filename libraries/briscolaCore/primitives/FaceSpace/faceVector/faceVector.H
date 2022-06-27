#ifndef faceVector_H
#define faceVector_H

#include "vector.H"
#include "FaceSpace.H"
#include "contiguous.H"

namespace Foam
{

typedef FaceSpace<vector> faceVector;

template<>
inline bool contiguous<faceVector>() { return true; }

}

#endif