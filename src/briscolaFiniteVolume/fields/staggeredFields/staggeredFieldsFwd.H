#ifndef staggeredFieldsFwd_H
#define staggeredFieldsFwd_H

#include "meshFieldFwd.H"
#include "staggeredLevelsFwd.H"

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

typedef meshField<label,staggered> staggeredLabelField;
typedef meshField<scalar,staggered> staggeredScalarField;
typedef meshField<faceScalar,staggered> staggeredFaceScalarField;
typedef meshField<vertexScalar,staggered> staggeredVertexScalarField;
typedef meshField<vector,staggered> staggeredVectorField;
typedef meshField<faceVector,staggered> staggeredFaceVectorField;
typedef meshField<vertexVector,staggered> staggeredVertexVectorField;
typedef meshField<tensor,staggered> staggeredTensorField;
typedef meshField<symmTensor,staggered> staggeredSymmTensorField;
typedef meshField<diagStencil,staggered> staggeredDiagStencilField;
typedef meshField<stencil,staggered> staggeredStencilField;

}

}

}

#endif
