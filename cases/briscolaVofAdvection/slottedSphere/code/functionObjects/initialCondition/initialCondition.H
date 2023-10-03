#ifndef functionObjects_initialCondition_H
#define functionObjects_initialCondition_H

#include "briscolaFunctionObject.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

class initialCondition
:
    public briscolaFunctionObject
{

    scalar TotalVolume_;
    scalar BoundError_;

public:

    // Runtime type information
    TypeName("initialCondition");

    // Construct from time and dictionary
    initialCondition
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Disallow default bitwise copy construction
    initialCondition(const initialCondition&) = delete;

    // Destructor
    virtual ~initialCondition();

    // Read the dictionary. Called at construction.
    virtual bool read(const dictionary&);

    // Code executed at the end of each time step loop
    virtual bool execute();

    // Code executed at the end of each time step loop, after execute()
    virtual bool write() { return true; }

    // Code executed at the end of the simulation
    virtual bool end();

    // Disallow default bitwise assignment
    void operator=(const initialCondition&) = delete;
};

}

}

}

}

#endif