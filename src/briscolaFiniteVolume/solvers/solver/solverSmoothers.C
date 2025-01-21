#include "solvers.H"
#include "addToRunTimeSelectionTable.H"

#include "RBGS.H"
#include "LEXGS.H"
#include "JAC.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeSmootherType(RBGS,diagStencil,scalar,colocated);
makeSmootherType(RBGS,diagStencil,scalar,staggered);
makeSmootherType(RBGS,diagStencil,vector,colocated);
makeSmootherType(RBGS,diagStencil,vector,staggered);

makeSmootherType(RBGS,stencil,scalar,colocated);
makeSmootherType(RBGS,stencil,scalar,staggered);
makeSmootherType(RBGS,stencil,vector,colocated);
makeSmootherType(RBGS,stencil,vector,staggered);


makeSmootherType(LEXGS,diagStencil,scalar,colocated);
makeSmootherType(LEXGS,diagStencil,scalar,staggered);
makeSmootherType(LEXGS,diagStencil,vector,colocated);
makeSmootherType(LEXGS,diagStencil,vector,staggered);

makeSmootherType(LEXGS,stencil,scalar,colocated);
makeSmootherType(LEXGS,stencil,scalar,staggered);
makeSmootherType(LEXGS,stencil,vector,colocated);
makeSmootherType(LEXGS,stencil,vector,staggered);


makeSmootherType(JAC,diagStencil,scalar,colocated);
makeSmootherType(JAC,diagStencil,scalar,staggered);
makeSmootherType(JAC,diagStencil,vector,colocated);
makeSmootherType(JAC,diagStencil,vector,staggered);

makeSmootherType(JAC,stencil,scalar,colocated);
makeSmootherType(JAC,stencil,scalar,staggered);
makeSmootherType(JAC,stencil,vector,colocated);
makeSmootherType(JAC,stencil,vector,staggered);


}

}

}
