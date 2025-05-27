#ifndef functionObjects_TaylorBubble_H
#define functionObjects_TaylorBubble_H

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

class TaylorBubble
:
    public briscolaFunctionObject
{
public:

    // Runtime type information
    TypeName("TaylorBubble");

    // Construct from time and dictionary
    TaylorBubble
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Disallow default bitwise copy construction
    TaylorBubble(const TaylorBubble&) = delete;

    // Destructor
    virtual ~TaylorBubble();

    // Read the dictionary. Called at construction.
    virtual bool read(const dictionary&);

    // Code executed at the end of each time step loop
    virtual bool execute();

    // Code executed at the end of each time step loop, after execute()
    virtual bool write() { return true; }

    // Code executed at the end of the simulation
    virtual bool end() { return true; }

    // Disallow default bitwise assignment
    void operator=(const TaylorBubble&) = delete;
};

}

}

}

}

#endif
