#ifndef functionObjects_Hysing_H
#define functionObjects_Hysing_H

#include "briscolaFunctionObject.H"
#include "OFstream.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

class Hysing
:
    public briscolaFunctionObject
{
    // Output file stream
    autoPtr<OFstream> filePtr_;

public:

    // Runtime type information
    TypeName("Hysing");

    // Construct from time and dictionary
    Hysing
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Disallow default bitwise copy construction
    Hysing(const Hysing&) = delete;

    // Destructor
    virtual ~Hysing();

    // Read the dictionary. Called at construction.
    virtual bool read(const dictionary&);

    // Code executed at the end of each time step loop
    virtual bool execute();

    // Code executed at the end of each time step loop, after execute()
    virtual bool write() { return true; }

    // Code executed at the end of the simulation
    virtual bool end() { return true; }

    // Disallow default bitwise assignment
    void operator=(const Hysing&) = delete;
};

}

}

}

}

#endif