#ifndef functionObjects_example_H
#define functionObjects_example_H

#include "briscolaFunctionObject.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

class example
:
    public briscolaFunctionObject
{
public:

    // Runtime type information
    TypeName("example");

    // Construct from time and dictionary
    example
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Disallow default bitwise copy construction
    example(const example&) = delete;

    // Destructor
    virtual ~example();

    // Read the dictionary. Called at construction.
    virtual bool read(const dictionary&);

    // Code executed at the end of each time step loop
    virtual bool execute();

    // Code executed at the end of each time step loop, after execute()
    virtual bool write();

    // Code executed at the end of the simulation
    virtual bool end();

    //- Return the list of fields required
    virtual wordList fields() const;

    // Disallow default bitwise assignment
    void operator=(const example&) = delete;
};

}

}

}

}

#endif
