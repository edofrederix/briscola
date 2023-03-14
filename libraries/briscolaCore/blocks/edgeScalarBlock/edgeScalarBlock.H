#ifndef edgeScalarBlock_H
#define edgeScalarBlock_H

#include "block.H"
#include "edgeScalar.H"

#define TEMPLATE
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

typedef block<edgeScalar> edgeScalarBlock;

BINARY_TYPE_OPERATOR(edgeScalar, scalar, edgeScalar, +, add)
BINARY_TYPE_OPERATOR(edgeScalar, scalar, edgeScalar, -, subtract)

BINARY_OPERATOR(edgeScalar, edgeScalar, edgeScalar, *, multiply)
BINARY_OPERATOR(edgeScalar, edgeScalar, edgeScalar, /, divide)

BINARY_TYPE_OPERATOR_SF(edgeScalar, edgeScalar, edgeScalar, /, divide)

}

}

#include "undefBlockFunctionsM.H"

#endif
