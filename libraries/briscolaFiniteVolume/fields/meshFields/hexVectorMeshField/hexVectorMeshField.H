#ifndef hexVectorMeshField_H
#define hexVectorMeshField_H

#include "hexVector.H"
#include "scalarMeshField.H"
#include "hexScalarMeshField.H"
#include "vectorMeshField.H"
#include "tensorMeshField.H"
#include "sphericalTensorMeshField.H"
#include "symmTensorMeshField.H"

#define TEMPLATE template<class MeshType>
#include "meshFieldFunctionsM.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

UNARY_FUNCTION(hexScalar, hexVector, magSqr)
UNARY_FUNCTION(hexScalar, hexVector, mag)
UNARY_FUNCTION(hexScalar, hexVector, cmptMax)
UNARY_FUNCTION(hexScalar, hexVector, cmptMin)
UNARY_FUNCTION(hexScalar, hexVector, cmptSum)
UNARY_FUNCTION(hexScalar, hexVector, cmptAv)
UNARY_FUNCTION(hexScalar, hexVector, cmptProduct)
UNARY_FUNCTION(hexVector, hexVector, cmptSqr)
UNARY_FUNCTION(hexVector, hexVector, cmptMag)

BINARY_FUNCTION(hexVector, hexVector, hexVector, cmptMultiply)
BINARY_FUNCTION(hexVector, hexVector, hexVector, cmptPow)
BINARY_FUNCTION(hexVector, hexVector, hexVector, cmptDivide)
BINARY_FUNCTION(hexVector, hexVector, hexVector, max)
BINARY_FUNCTION(hexVector, hexVector, hexVector, min)
BINARY_FUNCTION(hexVector, hexVector, hexVector, minMod)
BINARY_FUNCTION(hexVector, scalar, hexVector, dot)
BINARY_FUNCTION(hexVector, hexVector, scalar, dot)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "hexVectorMeshField.C"
#endif

#endif
