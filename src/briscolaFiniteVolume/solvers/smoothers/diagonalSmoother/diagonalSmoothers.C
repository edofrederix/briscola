#include "solvers.H"
#include "diagonalSmoother.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeSmootherType(diagonal,diagStencil,scalar,colocated);
makeSmootherType(diagonal,diagStencil,scalar,staggered);
makeSmootherType(diagonal,diagStencil,vector,colocated);
makeSmootherType(diagonal,diagStencil,vector,staggered);

makeSmootherType(diagonal,stencil,scalar,colocated);
makeSmootherType(diagonal,stencil,scalar,staggered);
makeSmootherType(diagonal,stencil,vector,colocated);
makeSmootherType(diagonal,stencil,vector,staggered);

}

}

}
