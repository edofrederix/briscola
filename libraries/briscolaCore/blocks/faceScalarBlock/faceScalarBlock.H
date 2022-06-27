#ifndef faceScalarBlock_H
#define faceScalarBlock_H

#include "block.H"
#include "faceScalar.H"

#define TEMPLATE
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

typedef block<faceScalar> faceScalarBlock;

BINARY_TYPE_OPERATOR(faceScalar, scalar, faceScalar, +, add)
BINARY_TYPE_OPERATOR(faceScalar, scalar, faceScalar, -, subtract)

BINARY_OPERATOR(faceScalar, faceScalar, faceScalar, *, multiply)
BINARY_OPERATOR(faceScalar, faceScalar, faceScalar, /, divide)

BINARY_TYPE_OPERATOR_SF(faceScalar, faceScalar, faceScalar, /, divide)

}

}

#include "undefBlockFunctionsM.H"

#endif
