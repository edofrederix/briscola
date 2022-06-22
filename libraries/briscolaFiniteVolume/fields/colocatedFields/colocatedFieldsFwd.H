#ifndef colocatedFieldsFwd_H
#define colocatedFieldsFwd_H

#include "colocatedLevelsFwd.H"

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

class colocated;

template<class Type, class MeshType>
class meshField;

typedef meshField<label,colocated> colocatedLabelField;
typedef meshField<scalar,colocated> colocatedScalarField;
typedef meshField<hexScalar,colocated> colocatedHexScalarField;
typedef meshField<vector,colocated> colocatedVectorField;
typedef meshField<hexVector,colocated> colocatedHexVectorField;
typedef meshField<tensor,colocated> colocatedTensorField;
typedef meshField<sphericalTensor,colocated> colocatedSphericalTensorField;
typedef meshField<symmTensor,colocated> colocatedSymmTensorField;
typedef meshField<diagTensor,colocated> colocatedDiagTensorField;
typedef meshField<symmStencil,colocated> colocatedSymmStencilField;
typedef meshField<diagStencil,colocated> colocatedDiagStencilField;
typedef meshField<stencil,colocated> colocatedStencilField;

}

}

}

#endif
