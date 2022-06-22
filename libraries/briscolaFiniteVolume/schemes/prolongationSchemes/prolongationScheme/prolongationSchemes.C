#include "prolongationSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "midPointProlongationScheme.H"
#include "linearProlongationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeProlongationScheme(scalar,colocated);
makeProlongationScheme(scalar,staggered);
makeProlongationScheme(vector,colocated);
makeProlongationScheme(vector,staggered);
makeProlongationScheme(tensor,colocated);
makeProlongationScheme(tensor,staggered);

makeProlongationSchemeType(midPoint,scalar,colocated);
makeProlongationSchemeType(midPoint,scalar,staggered);
makeProlongationSchemeType(midPoint,vector,colocated);
makeProlongationSchemeType(midPoint,vector,staggered);
makeProlongationSchemeType(midPoint,tensor,colocated);
makeProlongationSchemeType(midPoint,tensor,staggered);

makeProlongationSchemeType(linear,scalar,colocated);
makeProlongationSchemeType(linear,scalar,staggered);
makeProlongationSchemeType(linear,vector,colocated);
makeProlongationSchemeType(linear,vector,staggered);
makeProlongationSchemeType(linear,tensor,colocated);
makeProlongationSchemeType(linear,tensor,staggered);
}

}

}
