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

makeSmootherType(RBGS,stencil,scalar,colocated);
makeSmootherType(RBGS,stencil,scalar,staggered);
makeSmootherType(RBGS,stencil,vector,colocated);
makeSmootherType(RBGS,stencil,vector,staggered);

makeSmootherType(LEXGS,stencil,scalar,colocated);
makeSmootherType(LEXGS,stencil,scalar,staggered);
makeSmootherType(LEXGS,stencil,vector,colocated);
makeSmootherType(LEXGS,stencil,vector,staggered);

makeSmootherType(JAC,stencil,scalar,colocated);
makeSmootherType(JAC,stencil,scalar,staggered);
makeSmootherType(JAC,stencil,vector,colocated);
makeSmootherType(JAC,stencil,vector,staggered);


}

}

}
