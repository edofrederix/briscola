#ifndef lowerFaceScalarBlock_H
#define lowerFaceScalarBlock_H

#include "block.H"
#include "lowerFaceScalar.H"

#define TEMPLATE
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

typedef block<lowerFaceScalar> lowerFaceScalarBlock;

BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, scalar, +, add)
BINARY_TYPE_OPERATOR(lowerFaceScalar, scalar, lowerFaceScalar, +, add)

BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, scalar, -, subtract)
BINARY_TYPE_OPERATOR(lowerFaceScalar, scalar, lowerFaceScalar, -, subtract)

BINARY_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, /, divide)
BINARY_OPERATOR(lowerFaceScalar, faceScalar, lowerFaceScalar, /, divide)
BINARY_OPERATOR(lowerFaceScalar, lowerFaceScalar, faceScalar, /, divide)
BINARY_TYPE_OPERATOR(lowerFaceScalar, scalar, lowerFaceScalar, /, divide)

}

}

#include "undefBlockFunctionsM.H"

#endif
