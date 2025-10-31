#include "solvers.H"
#include "rbgsSmoother.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeSmootherType(rbgs,diagStencil,scalar,colocated);
makeSmootherType(rbgs,diagStencil,scalar,staggered);
makeSmootherType(rbgs,diagStencil,vector,colocated);
makeSmootherType(rbgs,diagStencil,vector,staggered);

makeSmootherType(rbgs,stencil,scalar,colocated);
makeSmootherType(rbgs,stencil,scalar,staggered);
makeSmootherType(rbgs,stencil,vector,colocated);
makeSmootherType(rbgs,stencil,vector,staggered);

}

}

}
