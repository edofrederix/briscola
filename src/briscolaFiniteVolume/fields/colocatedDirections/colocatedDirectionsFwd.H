#ifndef colocatedDirectionsFwd_H
#define colocatedDirectionsFwd_H

#include "meshDirectionFwd.H"

#include "label.H"
#include "scalar.H"
#include "faceScalar.H"
#include "vertexScalar.H"
#include "vector.H"
#include "faceVector.H"
#include "vertexVector.H"
#include "tensor.H"
#include "symmTensor.H"
#include "diagStencil.H"
#include "stencil.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class colocated;

typedef meshDirection<label,colocated> colocatedLabelDirection;
typedef meshDirection<scalar,colocated> colocatedScalarDirection;
typedef meshDirection<faceScalar,colocated> colocatedFaceScalarDirection;
typedef meshDirection<vertexScalar,colocated> colocatedVertexScalarDirection;
typedef meshDirection<vector,colocated> colocatedVectorDirection;
typedef meshDirection<faceVector,colocated> colocatedFaceVectorDirection;
typedef meshDirection<vertexVector,colocated> colocatedVertexVectorDirection;
typedef meshDirection<tensor,colocated> colocatedTensorDirection;
typedef meshDirection<symmTensor,colocated> colocatedSymmTensorDirection;
typedef meshDirection<diagStencil,colocated> colocatedDiagStencilDirection;
typedef meshDirection<stencil,colocated> colocatedStencilDirection;

}

}

}

#endif
