#ifndef geometryFun_H
#define geometryFun_H

#include "vector.H"

namespace Foam
{

namespace briscola
{

inline scalar tetVolume
(
    const vector& a,
    const vector& b,
    const vector& c,
    const vector& d
)
{
    return ((a-d) & ((c-d)^(b-d)))/6.0;
}

inline vector triAreaVector
(
    const vector& a,
    const vector& b,
    const vector& c
)
{
    return 0.5*((b-a) ^ (c-a));
}

inline scalar triArea
(
    const vector& a,
    const vector& b,
    const vector& c
)
{
    return Foam::mag(triAreaVector(a,b,c));
}

}

}

#endif
