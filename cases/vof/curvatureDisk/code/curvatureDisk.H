#ifndef functionObjects_curvatureDisk_H
#define functionObjects_curvatureDisk_H

#include "briscolaFunctionObject.H"
#include "meshFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

class curvatureDisk
:
    public briscolaFunctionObject
{
    const scalar R_;

public:

    // Runtime type information
    TypeName("curvatureDisk");

    // Construct from time and dictionary
    curvatureDisk
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    )
    :
        briscolaFunctionObject(name, runTime, dict),
        R_(0.5)
    {
        read(dict);
    }

    // Read the dictionary. Called at construction.
    virtual bool read(const dictionary&);

    // Code executed at the end of each time step loop
    virtual bool execute() { return true; }

    // Code executed at the end of each time step loop, after execute()
    virtual bool write() { return true; }

    // Code executed at the end of the simulation
    virtual bool end();
};

}

}

}

}

#endif