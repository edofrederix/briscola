#ifndef colocatedFieldsFwd_H
#define colocatedFieldsFwd_H

#include "meshFieldFwd.H"
#include "colocatedLevelsFwd.H"

#include "label.H"
#include "scalar.H"
#include "faceScalar.H"
#include "vertexScalar.H"
#include "vector.H"
#include "faceVector.H"
#include "vertexVector.H"
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

typedef meshField<label,colocated> colocatedLabelField;
typedef meshField<scalar,colocated> colocatedScalarField;
typedef meshField<faceScalar,colocated> colocatedFaceScalarField;
typedef meshField<vertexScalar,colocated> colocatedVertexScalarField;
typedef meshField<vector,colocated> colocatedVectorField;
typedef meshField<faceVector,colocated> colocatedFaceVectorField;
typedef meshField<vertexVector,colocated> colocatedVertexVectorField;
typedef meshField<tensor,colocated> colocatedTensorField;
typedef meshField<sphericalTensor,colocated> colocatedSphericalTensorField;
typedef meshField<symmTensor,colocated> colocatedSymmTensorField;
typedef meshField<diagTensor,colocated> colocatedDiagTensorField;
typedef meshField<diagStencil,colocated> colocatedDiagStencilField;
typedef meshField<stencil,colocated> colocatedStencilField;

}

}

}

#endif
