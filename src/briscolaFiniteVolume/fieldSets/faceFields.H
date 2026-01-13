#ifndef faceFields_H
#define faceFields_H

#include "faceFieldFwd.H"
#include "meshFieldSet.H"
#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

typedef faceField<scalar,colocated> colocatedScalarFaceField;
typedef faceField<vector,colocated> colocatedVectorFaceField;

typedef faceField<scalar,staggered> staggeredScalarFaceField;
typedef faceField<vector,staggered> staggeredVectorFaceField;

}

}

}

#endif
