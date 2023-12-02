#ifndef colocatedFieldsFwd_H
#define colocatedFieldsFwd_H

#include "colocatedLevelsFwd.H"

#include "label.H"
#include "scalar.H"
#include "faceScalar.H"
#include "lowerFaceScalar.H"
#include "edgeScalar.H"
#include "vertexScalar.H"
#include "vector.H"
#include "faceVector.H"
#include "edgeVector.H"
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

template<class Type, class MeshType>
class meshField;

typedef meshField<label,colocated> colocatedLabelField;
typedef meshField<scalar,colocated> colocatedScalarField;
typedef meshField<faceScalar,colocated> colocatedFaceScalarField;
typedef meshField<lowerFaceScalar,colocated> colocatedLowerFaceScalarField;
typedef meshField<edgeScalar,colocated> colocatedEdgeScalarField;
typedef meshField<vertexScalar,colocated> colocatedVertexScalarField;
typedef meshField<vector,colocated> colocatedVectorField;
typedef meshField<faceVector,colocated> colocatedFaceVectorField;
typedef meshField<edgeVector,colocated> colocatedEdgeVectorField;
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
