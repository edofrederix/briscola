#ifndef vertexLabelTensor_H
#define vertexLabelTensor_H

#include "labelTensor.H"
#include "VertexSpace.H"
#include "contiguous.H"

namespace Foam
{

typedef VertexSpace<labelTensor> vertexLabelTensor;

template<>
inline bool contiguous<vertexLabelTensor>() { return true; }

}

#endif