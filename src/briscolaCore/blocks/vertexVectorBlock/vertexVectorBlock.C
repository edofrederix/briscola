#include "vertexVectorBlock.H"

#define TEMPLATE template<int P>
#include "blockFunctionsM.C"

namespace Foam
{

namespace briscola
{

BINARY_TYPE_OPERATOR(vertexVector, vector, vertexVector, +, add)
BINARY_TYPE_OPERATOR(vertexVector, vector, vertexVector, -, subtract)

}

}

#include "undefBlockFunctionsM.H"