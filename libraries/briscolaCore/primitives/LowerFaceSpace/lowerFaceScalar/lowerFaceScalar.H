#ifndef lowerFaceScalar_H
#define lowerFaceScalar_H

#include "scalar.H"
#include "LowerFaceSpace.H"
#include "contiguous.H"

namespace Foam
{

typedef LowerFaceSpace<scalar> lowerFaceScalar;

template<>
inline bool contiguous<lowerFaceScalar>() { return true; }

}

#include "lowerFaceScalarI.H"

#endif
