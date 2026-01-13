#ifndef diagStencilFields_H
#define diagStencilFields_H

#include "diagStencilFieldFwd.H"
#include "meshFieldSet.H"
#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

typedef diagStencilField<scalar,colocated> colocatedScalarDiagStencilField;
typedef diagStencilField<vector,colocated> colocatedVectorDiagStencilField;

typedef diagStencilField<scalar,staggered> staggeredScalarDiagStencilField;
typedef diagStencilField<vector,staggered> staggeredVectorDiagStencilField;

}

}

}

#endif
