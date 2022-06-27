#ifndef colocatedLevelsFwd_H
#define colocatedLevelsFwd_H

#include "colocatedDirectionsFwd.H"

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

class colocated;

template<class Type, class MeshType>
class meshLevel;

typedef meshLevel<label,colocated> colocatedLabelLevel;
typedef meshLevel<scalar,colocated> colocatedScalarLevel;
typedef meshLevel<faceScalar,colocated> colocatedFaceScalarLevel;
typedef meshLevel<vector,colocated> colocatedVectorLevel;
typedef meshLevel<faceVector,colocated> colocatedFaceVectorLevel;
typedef meshLevel<tensor,colocated> colocatedTensorLevel;
typedef meshLevel<sphericalTensor,colocated> colocatedSphericalTensorLevel;
typedef meshLevel<symmTensor,colocated> colocatedSymmTensorLevel;
typedef meshLevel<diagTensor,colocated> colocatedDiagTensorLevel;
typedef meshLevel<diagStencil,colocated> colocatedDiagStencilLevel;
typedef meshLevel<stencil,colocated> colocatedStencilLevel;

}

}

}

#endif
