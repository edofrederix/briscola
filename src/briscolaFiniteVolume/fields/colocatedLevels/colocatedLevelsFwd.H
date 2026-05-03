#ifndef colocatedLevelsFwd_H
#define colocatedLevelsFwd_H

#include "meshLevelFwd.H"
#include "colocatedDirectionsFwd.H"

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

typedef meshLevel<label,colocated> colocatedLabelLevel;
typedef meshLevel<scalar,colocated> colocatedScalarLevel;
typedef meshLevel<faceScalar,colocated> colocatedFaceScalarLevel;
typedef meshLevel<vertexScalar,colocated> colocatedVertexScalarLevel;
typedef meshLevel<vector,colocated> colocatedVectorLevel;
typedef meshLevel<faceVector,colocated> colocatedFaceVectorLevel;
typedef meshLevel<vertexVector,colocated> colocatedVertexVectorLevel;
typedef meshLevel<tensor,colocated> colocatedTensorLevel;
typedef meshLevel<symmTensor,colocated> colocatedSymmTensorLevel;
typedef meshLevel<diagStencil,colocated> colocatedDiagStencilLevel;
typedef meshLevel<stencil,colocated> colocatedStencilLevel;

}

}

}

#endif
