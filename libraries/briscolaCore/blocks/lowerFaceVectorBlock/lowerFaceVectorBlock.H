#ifndef lowerFaceVectorBlock_H
#define lowerFaceVectorBlock_H

#include "lowerFaceScalarBlock.H"
#include "lowerFaceVector.H"
#include "vector.H"

#define TEMPLATE
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

typedef block<lowerFaceVector> lowerFaceVectorBlock;

BINARY_OPERATOR(lowerFaceVector, lowerFaceVector, lowerFaceScalar, /, divide)
BINARY_OPERATOR(lowerFaceVector, faceVector, lowerFaceScalar, /, divide)
BINARY_OPERATOR(lowerFaceVector, lowerFaceVector, faceScalar, /, divide)

UNARY_FUNCTION(lowerFaceScalar, lowerFaceVector, magSqr)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceVector, mag)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceVector, cmptMax)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceVector, cmptMin)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceVector, cmptSum)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceVector, cmptAv)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceVector, cmptProduct)
UNARY_FUNCTION(lowerFaceVector, lowerFaceVector, cmptSqr)
UNARY_FUNCTION(lowerFaceVector, lowerFaceVector, cmptMag)

BINARY_FUNCTION(lowerFaceVector, lowerFaceVector, lowerFaceVector, cmptMultiply)
BINARY_FUNCTION(lowerFaceVector, lowerFaceVector, lowerFaceVector, cmptPow)
BINARY_FUNCTION(lowerFaceVector, lowerFaceVector, lowerFaceVector, cmptDivide)
BINARY_FUNCTION(lowerFaceVector, lowerFaceVector, lowerFaceVector, max)
BINARY_FUNCTION(lowerFaceVector, lowerFaceVector, lowerFaceVector, min)
BINARY_FUNCTION(lowerFaceVector, lowerFaceVector, lowerFaceVector, minMod)

}

}

#include "undefBlockFunctionsM.H"

#endif
