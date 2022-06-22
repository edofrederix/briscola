#ifndef hexVectorMeshDirection_H
#define hexVectorMeshDirection_H

#include "hexVector.H"
#include "meshDirection.H"
#include "scalarMeshDirection.H"
#include "vectorMeshDirection.H"
#include "tensorMeshDirection.H"
#include "sphericalTensorMeshDirection.H"
#include "symmTensorMeshDirection.H"

#define TEMPLATE template<class MeshType>
#include "meshDirectionFunctionsM.H"

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

}

}

}

#include "undefBlockFunctionsM.H"


#ifdef NoRepository
    #include "hexVectorMeshDirection.C"
#endif

#endif
