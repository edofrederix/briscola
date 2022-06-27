#ifndef staggeredLevelsFwd_H
#define staggeredLevelsFwd_H

#include "staggeredDirectionsFwd.H"

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
class meshLevel;

typedef meshLevel<label,staggered> staggeredLabelLevel;
typedef meshLevel<scalar,staggered> staggeredScalarLevel;
typedef meshLevel<faceScalar,staggered> staggeredFaceScalarLevel;
typedef meshLevel<vector,staggered> staggeredVectorLevel;
typedef meshLevel<faceVector,staggered> staggeredFaceVectorLevel;
typedef meshLevel<tensor,staggered> staggeredTensorLevel;
typedef meshLevel<sphericalTensor,staggered> staggeredSphericalTensorLevel;
typedef meshLevel<symmTensor,staggered> staggeredSymmTensorLevel;
typedef meshLevel<diagTensor,staggered> staggeredDiagTensorLevel;
typedef meshLevel<diagStencil,staggered> staggeredDiagStencilLevel;
typedef meshLevel<stencil,staggered> staggeredStencilLevel;

}

}

}

#endif
