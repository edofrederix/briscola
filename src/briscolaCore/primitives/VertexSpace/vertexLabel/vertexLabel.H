#ifndef vertexLabel_H
#define vertexLabel_H

#include "label.H"
#include "VertexSpace.H"
#include "contiguous.H"

namespace Foam
{

typedef VertexSpace<label> vertexLabel;

template<>
inline bool contiguous<vertexLabel>() { return true; }

}

#endif