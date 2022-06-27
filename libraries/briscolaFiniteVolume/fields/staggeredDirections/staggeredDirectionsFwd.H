#ifndef staggeredDirectionsFwd_H
#define staggeredDirectionsFwd_H

#include "label.H"
#include "scalar.H"
#include "faceScalar.H"
#include "vector.H"
#include "faceVector.H"
#include "tensor.H"
#include "sphericalTensor.H"
#include "symmTensor.H"
#include "diagTensor.H"
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
typedef meshDirection<faceScalar,staggered> staggeredFaceScalarDirection;
typedef meshDirection<vector,staggered> staggeredVectorDirection;
typedef meshDirection<faceVector,staggered> staggeredFaceVectorDirection;
typedef meshDirection<tensor,staggered> staggeredTensorDirection;
typedef meshDirection<sphericalTensor,staggered> staggeredSphericalTensorDirection;
typedef meshDirection<symmTensor,staggered> staggeredSymmTensorDirection;
typedef meshDirection<diagTensor,staggered> staggeredDiagTensorDirection;
typedef meshDirection<diagStencil,staggered> staggeredDiagStencilDirection;
typedef meshDirection<stencil,staggered> staggeredStencilDirection;

}

}

}

#endif
