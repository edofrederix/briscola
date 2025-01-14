#include "ddtSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "EulerDdtScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeDdtScheme(scalar,colocated)
makeDdtScheme(scalar,staggered)

makeDdtScheme(vector,colocated)
makeDdtScheme(vector,staggered)

makeDdtSchemeType(Euler,scalar,colocated)
makeDdtSchemeType(Euler,scalar,staggered)

makeDdtSchemeType(Euler,vector,colocated)
makeDdtSchemeType(Euler,vector,staggered)

}

}

}
