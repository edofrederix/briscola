#ifndef vectorBlock_H
#define vectorBlock_H

#include "scalarBlock.H"
#include "vector.H"

namespace Foam
{

namespace briscola
{

typedef block<vector> vectorBlock;

}

// Resolve conversion ambiguities

inline vector operator/(const vector& v, scalar s)
{
    return vector(v.x()/s, v.y()/s, v.z()/s);
}

}

#endif
