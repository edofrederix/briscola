#ifndef faceLabelTensor_H
#define faceLabelTensor_H

#include "labelTensor.H"
#include "FaceSpace.H"
#include "contiguous.H"

namespace Foam
{

typedef FaceSpace<labelTensor> faceLabelTensor;

template<>
inline bool contiguous<faceLabelTensor>() { return true; }

}

#endif
