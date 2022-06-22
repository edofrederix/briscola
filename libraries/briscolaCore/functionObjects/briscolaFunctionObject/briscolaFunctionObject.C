#include "briscolaFunctionObject.H"
#include "Time.H"

namespace Foam
{

namespace briscola
{

briscolaFunctionObject::briscolaFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    runTime_(runTime),
    dict_(dict)
{}

briscolaFunctionObject::~briscolaFunctionObject()
{}

}

}