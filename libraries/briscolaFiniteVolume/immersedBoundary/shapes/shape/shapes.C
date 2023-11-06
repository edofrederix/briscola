#include "shapes.H"
#include "addToRunTimeSelectionTable.H"

#include "cylinder.H"
#include "sphere.H"

// Also compile non-templated shape
#include "shape.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(shape, 0);
defineRunTimeSelectionTable(shape, dictionary);

makeShapeType(cylinder)
makeShapeType(sphere)

}

}

}