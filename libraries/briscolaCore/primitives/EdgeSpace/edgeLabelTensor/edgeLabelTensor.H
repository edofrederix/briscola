#ifndef edgeLabelTensor_H
#define edgeLabelTensor_H

#include "labelTensor.H"
#include "EdgeSpace.H"
#include "contiguous.H"

namespace Foam
{

typedef EdgeSpace<labelTensor> edgeLabelTensor;

template<>
inline bool contiguous<edgeLabelTensor>() { return true; }

}

#endif
