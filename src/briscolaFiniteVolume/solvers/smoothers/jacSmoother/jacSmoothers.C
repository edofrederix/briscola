#include "solvers.H"
#include "jacSmoother.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeSmootherType(jac,diagStencil,scalar,colocated);
makeSmootherType(jac,diagStencil,scalar,staggered);
makeSmootherType(jac,diagStencil,vector,colocated);
makeSmootherType(jac,diagStencil,vector,staggered);

makeSmootherType(jac,stencil,scalar,colocated);
makeSmootherType(jac,stencil,scalar,staggered);
makeSmootherType(jac,stencil,vector,colocated);
makeSmootherType(jac,stencil,vector,staggered);

}

}

}
