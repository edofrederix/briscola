#ifndef colocatedDirectionsFwd_H
#define colocatedDirectionsFwd_H

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
class meshDirection;

typedef meshDirection<label,colocated> colocatedLabelDirection;
typedef meshDirection<scalar,colocated> colocatedScalarDirection;
typedef meshDirection<hexScalar,colocated> colocatedHexScalarDirection;
typedef meshDirection<vector,colocated> colocatedVectorDirection;
typedef meshDirection<hexVector,colocated> colocatedHexVectorDirection;
typedef meshDirection<tensor,colocated> colocatedTensorDirection;
typedef meshDirection<sphericalTensor,colocated> colocatedSphericalTensorDirection;
typedef meshDirection<symmTensor,colocated> colocatedSymmTensorDirection;
typedef meshDirection<diagTensor,colocated> colocatedDiagTensorDirection;
typedef meshDirection<symmStencil,colocated> colocatedSymmStencilDirection;
typedef meshDirection<diagStencil,colocated> colocatedDiagStencilDirection;
typedef meshDirection<stencil,colocated> colocatedStencilDirection;

}

}

}

#endif
