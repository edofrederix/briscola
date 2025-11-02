#ifndef functionObjects_flowOverCylinder_H
#define functionObjects_flowOverCylinder_H

#include "briscolaFunctionObject.H"
#include "pointInterpolator.H"
#include "OFstream.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

class flowOverCylinder
:
    public briscolaFunctionObject
{
    // Output file stream
    autoPtr<OFstream> filePtr_;

    // Interpolator
    autoPtr<pointInterpolator<colocated>> interpPtr_;

public:

    // Runtime type information
    TypeName("flowOverCylinder");

    // Construct from time and dictionary
    flowOverCylinder
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Read the dictionary. Called at construction.
    virtual bool read(const dictionary&) { return true; }

    // Code executed at the end of each time step loop
    virtual bool execute();

    // Code executed at the end of each time step loop, after execute()
    virtual bool write() { return true; }

    // Code executed at the end of the simulation
    virtual bool end() { return true; }
};

}

}

}

}

#endif