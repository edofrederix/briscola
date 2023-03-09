#ifndef vertexScalarBlock_H
#define vertexScalarBlock_H

#include "block.H"
#include "vertexScalar.H"

#define TEMPLATE
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

typedef block<vertexScalar> vertexScalarBlock;

BINARY_TYPE_OPERATOR(vertexScalar, scalar, vertexScalar, +, add)
BINARY_TYPE_OPERATOR(vertexScalar, scalar, vertexScalar, -, subtract)

BINARY_OPERATOR(vertexScalar, vertexScalar, vertexScalar, *, multiply)
BINARY_OPERATOR(vertexScalar, vertexScalar, vertexScalar, /, divide)

BINARY_TYPE_OPERATOR_SF(vertexScalar, vertexScalar, vertexScalar, /, divide)

}

}

#include "undefBlockFunctionsM.H"

#endif
