#ifndef stencilFields_H
#define stencilFields_H

#include "stencilFieldFwd.H"
#include "meshFieldSet.H"
#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

typedef stencilField<scalar,colocated> colocatedScalarStencilField;
typedef stencilField<vector,colocated> colocatedVectorStencilField;

typedef stencilField<scalar,staggered> staggeredScalarStencilField;
typedef stencilField<vector,staggered> staggeredVectorStencilField;

}

}

}

#endif
