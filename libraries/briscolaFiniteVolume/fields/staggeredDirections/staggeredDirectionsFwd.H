#ifndef staggeredDirectionsFwd_H
#define staggeredDirectionsFwd_H

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
class meshDirection;

typedef meshDirection<label,staggered> staggeredLabelDirection;
typedef meshDirection<scalar,staggered> staggeredScalarDirection;
typedef meshDirection<hexScalar,staggered> staggeredHexScalarDirection;
typedef meshDirection<vector,staggered> staggeredVectorDirection;
typedef meshDirection<hexVector,staggered> staggeredHexVectorDirection;
typedef meshDirection<tensor,staggered> staggeredTensorDirection;
typedef meshDirection<sphericalTensor,staggered> staggeredSphericalTensorDirection;
typedef meshDirection<symmTensor,staggered> staggeredSymmTensorDirection;
typedef meshDirection<diagTensor,staggered> staggeredDiagTensorDirection;
typedef meshDirection<symmStencil,staggered> staggeredSymmStencilDirection;
typedef meshDirection<diagStencil,staggered> staggeredDiagStencilDirection;
typedef meshDirection<stencil,staggered> staggeredStencilDirection;

}

}

}

#endif
