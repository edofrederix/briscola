#ifndef colocatedLevelsFwd_H
#define colocatedLevelsFwd_H

#include "colocatedDirectionsFwd.H"

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
class meshLevel;

typedef meshLevel<label,colocated> colocatedLabelLevel;
typedef meshLevel<scalar,colocated> colocatedScalarLevel;
typedef meshLevel<hexScalar,colocated> colocatedHexScalarLevel;
typedef meshLevel<vector,colocated> colocatedVectorLevel;
typedef meshLevel<hexVector,colocated> colocatedHexVectorLevel;
typedef meshLevel<tensor,colocated> colocatedTensorLevel;
typedef meshLevel<sphericalTensor,colocated> colocatedSphericalTensorLevel;
typedef meshLevel<symmTensor,colocated> colocatedSymmTensorLevel;
typedef meshLevel<diagTensor,colocated> colocatedDiagTensorLevel;
typedef meshLevel<symmStencil,colocated> colocatedSymmStencilLevel;
typedef meshLevel<diagStencil,colocated> colocatedDiagStencilLevel;
typedef meshLevel<stencil,colocated> colocatedStencilLevel;

}

}

}

#endif
