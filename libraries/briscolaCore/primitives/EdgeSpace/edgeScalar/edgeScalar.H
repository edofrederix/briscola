#ifndef edgeScalar_H
#define edgeScalar_H

#include "scalar.H"
#include "EdgeSpace.H"
#include "contiguous.H"

namespace Foam
{

typedef EdgeSpace<scalar> edgeScalar;

template<>
inline bool contiguous<edgeScalar>() { return true; }

// Specialization of CellSpaceCast for stencils

template<>
class CellSpaceCast<edgeScalar, stencil, scalar>
{
public:
    typedef stencil type;
};

template<>
class CellSpaceCast<stencil, edgeScalar, scalar>
{
public:
    typedef stencil type;
};

}

#endif