#include "flowOverCylinder.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include <fstream>

#include "meshFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

defineTypeNameAndDebug(flowOverCylinder, 0);

addToRunTimeSelectionTable
(
    functionObject,
    flowOverCylinder,
    dictionary
);

flowOverCylinder::flowOverCylinder
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    briscolaFunctionObject(name, runTime, dict)
{
    read(dict);

    // Write data to data.txt

    if (Pstream::master())
    {
        filePtr_.reset(new OFstream("data.txt"));
    }
}

flowOverCylinder::~flowOverCylinder()
{}

bool flowOverCylinder::read(const dictionary& dict)
{
    return true;
}

bool flowOverCylinder::execute()
{
    // Find the velocity field

    const colocatedVectorField& U =
        runTime_.lookupObjectRef<colocatedVectorField>("U");

    if (!interpPtr_.valid())
    {
        // Setup the point interpolator if not yet done. Sample the velocity in
        // a point that's at x = 4D in the wake of the cylinder.

        vectorList points(1, vector(4,0,0));

        interpPtr_.reset
        (
            pointInterpolator<colocated>::New(U.fvMsh(), points, "linear").ptr()
        );
    }

    // Sample the list of points (just one)

    vectorList data(this->interpPtr_->operator()(U));

    // Let the master process write the outcome to the output file

    if (Pstream::master())
        filePtr_()
            << runTime_.time().value()
            << " " << data[0].x()
            << " " << data[0].y()
            << " " << data[0].z() << nl;

    return true;
}

}

}

}

}
