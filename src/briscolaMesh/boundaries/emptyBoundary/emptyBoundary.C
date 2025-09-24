#include "emptyBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(emptyBoundary, 0);
addToRunTimeSelectionTable(boundary, emptyBoundary, dictionary);

emptyBoundary::emptyBoundary(const mesh& msh, const dictionary& dict)
:
    boundary(msh, dict)
{
    // Empty patches are extended patches

    this->extend();
}

emptyBoundary::emptyBoundary(const emptyBoundary& pp)
:
    boundary(pp)
{}

}

}
