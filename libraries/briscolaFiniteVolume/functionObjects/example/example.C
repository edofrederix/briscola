#include "example.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

defineTypeNameAndDebug(example, 0);

addToRunTimeSelectionTable
(
    functionObject,
    example,
    dictionary
);

example::example
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    briscolaFunctionObject(name, runTime, dict)
{}

example::~example()
{}

bool example::execute()
{
    Info<< "Call to execute() at t = "
        << runTime_.time().value() << endl;

    return true;
}

bool example::write()
{
    Info<< "Call to write() at t = "
        << runTime_.time().value() << endl;

    return true;
}

bool example::end()
{
    Info<< "Call to end() at t = "
        << runTime_.time().value() << endl;

    return true;
}

}

}

}

}
