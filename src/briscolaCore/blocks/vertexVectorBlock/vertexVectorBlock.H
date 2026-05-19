#ifndef vertexVectorBlock_H
#define vertexVectorBlock_H

#include "vertexScalarBlock.H"
#include "vertexVector.H"
#include "vector.H"

#define TEMPLATE template<int P>
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

BINARY_TYPE_OPERATOR(vertexVector, vector, vertexVector, +, add)
BINARY_TYPE_OPERATOR(vertexVector, vector, vertexVector, -, subtract)

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "vertexVectorBlock.C"
#endif

#endif
