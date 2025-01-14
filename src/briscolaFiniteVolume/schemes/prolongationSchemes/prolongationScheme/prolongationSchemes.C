#include "prolongationSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "linearProlongationScheme.H"
#include "injectionProlongationScheme.H"

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

makeProlongationSchemeType(linear,scalar,colocated);
makeProlongationSchemeType(linear,scalar,staggered);
makeProlongationSchemeType(linear,vector,colocated);
makeProlongationSchemeType(linear,vector,staggered);
makeProlongationSchemeType(linear,tensor,colocated);
makeProlongationSchemeType(linear,tensor,staggered);

makeProlongationSchemeType(injection,scalar,colocated);
makeProlongationSchemeType(injection,scalar,staggered);
makeProlongationSchemeType(injection,vector,colocated);
makeProlongationSchemeType(injection,vector,staggered);
makeProlongationSchemeType(injection,tensor,colocated);
makeProlongationSchemeType(injection,tensor,staggered);

}

}

}
