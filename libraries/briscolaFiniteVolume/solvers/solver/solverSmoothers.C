#include "solvers.H"
#include "addToRunTimeSelectionTable.H"

#include "RBGS.H"
#include "LEXGS.H"
#include "JAC.H"

#include "symmRBGS.H"
#include "symmLEXGS.H"
#include "symmJAC.H"

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

makeSmootherType(RBGS,symmStencil,scalar,colocated);
makeSmootherType(RBGS,symmStencil,scalar,staggered);
makeSmootherType(RBGS,symmStencil,vector,colocated);
makeSmootherType(RBGS,symmStencil,vector,staggered);

makeSmootherType(symmRBGS,stencil,scalar,colocated);
makeSmootherType(symmRBGS,stencil,scalar,staggered);
makeSmootherType(symmRBGS,stencil,vector,colocated);
makeSmootherType(symmRBGS,stencil,vector,staggered);

makeSmootherType(symmRBGS,symmStencil,scalar,colocated);
makeSmootherType(symmRBGS,symmStencil,scalar,staggered);
makeSmootherType(symmRBGS,symmStencil,vector,colocated);
makeSmootherType(symmRBGS,symmStencil,vector,staggered);


makeSmootherType(LEXGS,stencil,scalar,colocated);
makeSmootherType(LEXGS,stencil,scalar,staggered);
makeSmootherType(LEXGS,stencil,vector,colocated);
makeSmootherType(LEXGS,stencil,vector,staggered);

makeSmootherType(LEXGS,symmStencil,scalar,colocated);
makeSmootherType(LEXGS,symmStencil,scalar,staggered);
makeSmootherType(LEXGS,symmStencil,vector,colocated);
makeSmootherType(LEXGS,symmStencil,vector,staggered);

makeSmootherType(symmLEXGS,stencil,scalar,colocated);
makeSmootherType(symmLEXGS,stencil,scalar,staggered);
makeSmootherType(symmLEXGS,stencil,vector,colocated);
makeSmootherType(symmLEXGS,stencil,vector,staggered);

makeSmootherType(symmLEXGS,symmStencil,scalar,colocated);
makeSmootherType(symmLEXGS,symmStencil,scalar,staggered);
makeSmootherType(symmLEXGS,symmStencil,vector,colocated);
makeSmootherType(symmLEXGS,symmStencil,vector,staggered);


makeSmootherType(JAC,stencil,scalar,colocated);
makeSmootherType(JAC,stencil,scalar,staggered);
makeSmootherType(JAC,stencil,vector,colocated);
makeSmootherType(JAC,stencil,vector,staggered);

makeSmootherType(JAC,symmStencil,scalar,colocated);
makeSmootherType(JAC,symmStencil,scalar,staggered);
makeSmootherType(JAC,symmStencil,vector,colocated);
makeSmootherType(JAC,symmStencil,vector,staggered);

makeSmootherType(symmJAC,stencil,scalar,colocated);
makeSmootherType(symmJAC,stencil,scalar,staggered);
makeSmootherType(symmJAC,stencil,vector,colocated);
makeSmootherType(symmJAC,stencil,vector,staggered);

makeSmootherType(symmJAC,symmStencil,scalar,colocated);
makeSmootherType(symmJAC,symmStencil,scalar,staggered);
makeSmootherType(symmJAC,symmStencil,vector,colocated);
makeSmootherType(symmJAC,symmStencil,vector,staggered);

}

}

}
