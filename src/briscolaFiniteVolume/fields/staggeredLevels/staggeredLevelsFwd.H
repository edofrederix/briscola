#ifndef staggeredLevelsFwd_H
#define staggeredLevelsFwd_H

#include "meshLevelFwd.H"
#include "staggeredDirectionsFwd.H"

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

class staggered;

typedef meshLevel<label,staggered> staggeredLabelLevel;
typedef meshLevel<scalar,staggered> staggeredScalarLevel;
typedef meshLevel<faceScalar,staggered> staggeredFaceScalarLevel;
typedef meshLevel<vertexScalar,staggered> staggeredVertexScalarLevel;
typedef meshLevel<vector,staggered> staggeredVectorLevel;
typedef meshLevel<faceVector,staggered> staggeredFaceVectorLevel;
typedef meshLevel<vertexVector,staggered> staggeredVertexVectorLevel;
typedef meshLevel<tensor,staggered> staggeredTensorLevel;
typedef meshLevel<symmTensor,staggered> staggeredSymmTensorLevel;
typedef meshLevel<diagStencil,staggered> staggeredDiagStencilLevel;
typedef meshLevel<stencil,staggered> staggeredStencilLevel;

}

}

}

#endif
