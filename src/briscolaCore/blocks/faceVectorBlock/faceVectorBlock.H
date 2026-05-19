#ifndef faceVectorBlock_H
#define faceVectorBlock_H

#include "faceScalarBlock.H"
#include "faceVector.H"
#include "vector.H"

#define TEMPLATE template<int P>
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

BINARY_OPERATOR(faceVector, faceVector, faceScalar, /, divide)

BINARY_TYPE_FUNCTION(faceVector, faceVector, vector, max)
BINARY_TYPE_FUNCTION(faceVector, vector, faceVector, max)

BINARY_TYPE_FUNCTION(faceVector, faceVector, vector, min)
BINARY_TYPE_FUNCTION(faceVector, vector, faceVector, min)

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "faceVectorBlock.C"
#endif

#endif
