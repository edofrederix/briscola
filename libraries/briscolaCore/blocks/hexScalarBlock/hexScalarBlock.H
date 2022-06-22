#ifndef hexScalarBlock_H
#define hexScalarBlock_H

#include "block.H"
#include "hexScalar.H"

#define TEMPLATE
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

typedef block<hexScalar> hexScalarBlock;

BINARY_TYPE_OPERATOR(hexScalar, hexScalar, hexScalar, +, add)
BINARY_TYPE_OPERATOR(hexScalar, hexScalar, hexScalar, -, subtract)

BINARY_OPERATOR(hexScalar, hexScalar, hexScalar, *, multiply)
BINARY_OPERATOR(hexScalar, hexScalar, hexScalar, /, divide)

BINARY_TYPE_OPERATOR_SF(hexScalar, hexScalar, hexScalar, /, divide)

}

}

#include "undefBlockFunctionsM.H"

#endif
