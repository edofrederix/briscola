#include "immersedBoundaryShapes.H"
#include "addToRunTimeSelectionTable.H"

#include "immersedBoundaryCylinder.H"
#include "immersedBoundarySphere.H"
#include "immersedBoundarySlab.H"

// Also compile non-templated shape
#include "immersedBoundaryShape.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(immersedBoundaryShape, 0);
defineRunTimeSelectionTable(immersedBoundaryShape, dictionary);

makeShapeType(immersedBoundaryCylinder)
makeShapeType(immersedBoundarySphere)
makeShapeType(immersedBoundarySlab)

}

}

}