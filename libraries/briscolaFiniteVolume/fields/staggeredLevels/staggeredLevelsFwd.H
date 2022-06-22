#ifndef staggeredLevelsFwd_H
#define staggeredLevelsFwd_H

#include "staggeredDirectionsFwd.H"

#include "label.H"
#include "scalar.H"
#include "hexScalar.H"
#include "vector.H"
#include "hexVector.H"
#include "tensor.H"
#include "sphericalTensor.H"
#include "symmTensor.H"
#include "diagTensor.H"
#include "symmStencil.H"
#include "diagStencil.H"
#include "stencil.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class staggered;

template<class Type, class MeshType>
class meshLevel;

typedef meshLevel<label,staggered> staggeredLabelLevel;
typedef meshLevel<scalar,staggered> staggeredScalarLevel;
typedef meshLevel<hexScalar,staggered> staggeredHexScalarLevel;
typedef meshLevel<vector,staggered> staggeredVectorLevel;
typedef meshLevel<hexVector,staggered> staggeredHexVectorLevel;
typedef meshLevel<tensor,staggered> staggeredTensorLevel;
typedef meshLevel<sphericalTensor,staggered> staggeredSphericalTensorLevel;
typedef meshLevel<symmTensor,staggered> staggeredSymmTensorLevel;
typedef meshLevel<diagTensor,staggered> staggeredDiagTensorLevel;
typedef meshLevel<symmStencil,staggered> staggeredSymmStencilLevel;
typedef meshLevel<diagStencil,staggered> staggeredDiagStencilLevel;
typedef meshLevel<stencil,staggered> staggeredStencilLevel;

}

}

}

#endif
