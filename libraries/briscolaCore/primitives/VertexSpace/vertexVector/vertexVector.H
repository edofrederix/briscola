#ifndef vertexVector_H
#define vertexVector_H

#include "vector.H"
#include "VertexSpace.H"
#include "contiguous.H"

namespace Foam
{

typedef VertexSpace<vector> vertexVector;

template<>
inline bool contiguous<vertexVector>() { return true; }

// Given a point in a general quadrilaterally-faced hexahedron defined by eight
// vertices, compute the interpolation weights to that point. Optionally, only
// calculate the local coordinates if the point is inside the hexahedron,
// allowing for early rejection returning -(1,1,1) if the point is not in the
// hexahedron.
vector interpolationWeights
(
    const vector& point,
    const vertexVector& vertices,
    const bool insideOnly = false
);

}

#endif