#ifndef vertexVectorBlock_H
#define vertexVectorBlock_H

#include "vertexScalarBlock.H"
#include "vertexVector.H"
#include "vector.H"

#define TEMPLATE
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

typedef block<vertexVector> vertexVectorBlock;

BINARY_TYPE_OPERATOR(vertexVector, vector, vertexVector, +, add)
BINARY_TYPE_OPERATOR(vertexVector, vector, vertexVector, -, subtract)

}

}

#include "undefBlockFunctionsM.H"

#endif
