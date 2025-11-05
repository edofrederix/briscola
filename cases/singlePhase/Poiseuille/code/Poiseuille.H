#ifndef functionObjects_Poiseuille_H
#define functionObjects_Poiseuille_H

#include "briscolaFunctionObject.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

class Poiseuille
:
    public briscolaFunctionObject
{
public:

    // Runtime type information
    TypeName("Poiseuille");

    // Construct from time and dictionary
    Poiseuille
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