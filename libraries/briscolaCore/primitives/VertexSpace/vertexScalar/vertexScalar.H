#ifndef vertexScalar_H
#define vertexScalar_H

#include "scalar.H"
#include "VertexSpace.H"
#include "contiguous.H"

namespace Foam
{

typedef VertexSpace<scalar> vertexScalar;

template<>
inline bool contiguous<vertexScalar>() { return true; }

// Specialization of CellSpaceCast for stencils

template<>
class CellSpaceCast<vertexScalar, stencil, scalar>
{
public:
    typedef stencil type;
};

template<>
class CellSpaceCast<stencil, vertexScalar, scalar>
{
public:
    typedef stencil type;
};

}

#endif