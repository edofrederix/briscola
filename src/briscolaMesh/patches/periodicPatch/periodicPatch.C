#include "periodicPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "IOobject.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(periodicPatch, 0);
addToRunTimeSelectionTable(patch, periodicPatch, dictionary);

periodicPatch::periodicPatch
(
    const geometry& g,
    const label num,
    const word name,
    const dictionary& dict
)
:
    patch(g, num, name, dict),
    neighbor_(dict.lookup("neighbor"))
{
    // If a sub-number other than zero is found, the periodic patch is invalid
    // because it may not consume more than one face of the same brick.

    const label subNum(atoi(IOobject::group(name).c_str()));

    if (subNum)
        FatalErrorInFunction
            << "A periodic patch cannot have two faces in the same brick"
            << endl << abort(FatalError);
}

periodicPatch::periodicPatch(const periodicPatch& p)
:
    patch(p)
{}

periodicPatch::~periodicPatch()
{}

}

}
