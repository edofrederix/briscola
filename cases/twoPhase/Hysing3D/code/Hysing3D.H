#ifndef functionObjects_Hysing3D_H
#define functionObjects_Hysing3D_H

#include "briscolaFunctionObject.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

class Hysing3D
:
    public briscolaFunctionObject
{
public:

    // Runtime type information
    TypeName("Hysing3D");

    // Construct from time and dictionary
    Hysing3D
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    )
    :
        briscolaFunctionObject(name, runTime, dict)
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
    virtual bool end() { return true; }
};

}

}

}

}

#endif
