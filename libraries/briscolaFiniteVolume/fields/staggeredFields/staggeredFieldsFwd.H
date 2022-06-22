#ifndef staggeredFieldsFwd_H
#define staggeredFieldsFwd_H

#include "staggeredLevelsFwd.H"

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
class meshField;

typedef meshField<label,staggered> staggeredLabelField;
typedef meshField<scalar,staggered> staggeredScalarField;
typedef meshField<hexScalar,staggered> staggeredHexScalarField;
typedef meshField<vector,staggered> staggeredVectorField;
typedef meshField<hexVector,staggered> staggeredHexVectorField;
typedef meshField<tensor,staggered> staggeredTensorField;
typedef meshField<sphericalTensor,staggered> staggeredSphericalTensorField;
typedef meshField<symmTensor,staggered> staggeredSymmTensorField;
typedef meshField<diagTensor,staggered> staggeredDiagTensorField;
typedef meshField<symmStencil,staggered> staggeredSymmStencilField;
typedef meshField<diagStencil,staggered> staggeredDiagStencilField;
typedef meshField<stencil,staggered> staggeredStencilField;

}

}

}

#endif
