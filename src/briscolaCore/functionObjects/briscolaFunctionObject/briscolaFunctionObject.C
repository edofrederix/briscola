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
    functionObject(name, runTime),
    runTime_(runTime),
    name_(name),
    dict_(dict)
{
    read(dict);
}

briscolaFunctionObject::~briscolaFunctionObject()
{}

bool briscolaFunctionObject::read(const dictionary& dict)
{
    return functionObject::read(dict);
}

}

}