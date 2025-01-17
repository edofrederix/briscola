#ifndef faceLabelVector_H
#define faceLabelVector_H

#include "labelVector.H"
#include "FaceSpace.H"
#include "contiguous.H"

namespace Foam
{

typedef FaceSpace<labelVector> faceLabelVector;

template<>
inline bool contiguous<faceLabelVector>() { return true; }

}

#endif
