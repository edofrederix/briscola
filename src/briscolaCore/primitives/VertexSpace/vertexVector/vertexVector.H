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

namespace briscola
{

// Volume of a hexahedron computed from a tetrahedron decomposition
inline scalar hexVolume(const vertexVector& v)
{
    return
        tetVolume(v[4], v[0], v[5], v[6])
      + tetVolume(v[1], v[5], v[0], v[3])
      + tetVolume(v[2], v[6], v[3], v[0])
      + tetVolume(v[7], v[3], v[6], v[5])
      + tetVolume(v[6], v[3], v[0], v[5]);
}

// Given a point in a general quadrilaterally-faced hexahedron defined by eight
// vertices, compute the interpolation weights to that point. Optionally, only
// calculate the local coordinates if the point is inside the hexahedron,
// allowing for early rejection returning -(1,1,1) if the point is not in the
// hexahedron. If convergence fails and fatal is switched off, the function
// returns -(1,1,1) without error.
vector interpolationWeights
(
    const vector& point,
    const vertexVector& vertices,
    const bool insideOnly = false,
    const bool fatal = true,
    const scalar tol = 1e-12
);

}

}

#endif