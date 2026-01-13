#include "solvers.H"
#include "lexgsSmoother.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeSmootherType(lexgs,diagStencil,scalar,colocated);
makeSmootherType(lexgs,diagStencil,scalar,staggered);
makeSmootherType(lexgs,diagStencil,vector,colocated);
makeSmootherType(lexgs,diagStencil,vector,staggered);

makeSmootherType(lexgs,stencil,scalar,colocated);
makeSmootherType(lexgs,stencil,scalar,staggered);
makeSmootherType(lexgs,stencil,vector,colocated);
makeSmootherType(lexgs,stencil,vector,staggered);

}

}

}
