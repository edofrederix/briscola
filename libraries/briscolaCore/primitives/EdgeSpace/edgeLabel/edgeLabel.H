#ifndef edgeLabel_H
#define edgeLabel_H

#include "label.H"
#include "EdgeSpace.H"
#include "contiguous.H"

namespace Foam
{

typedef EdgeSpace<label> edgeLabel;

template<>
inline bool contiguous<edgeLabel>() { return true; }

}

#endif