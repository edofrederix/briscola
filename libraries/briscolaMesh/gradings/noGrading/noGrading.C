#include "noGrading.H"
#include "geometry.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace gradings
{

defineTypeNameAndDebug(noGrading, 0);
addToRunTimeSelectionTable(grading, noGrading, dictionary);

noGrading::noGrading(const brick& b)
:
    grading(b)
{}

noGrading::noGrading(const noGrading& g)
:
    grading(g)
{}

}

}

}
