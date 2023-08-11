#ifndef briscolaFunctionObject_H
#define briscolaFunctionObject_H

#include "functionObject.H"

namespace Foam
{

namespace briscola
{

class briscolaFunctionObject
:
    public functionObject
{
protected:

    // Reference to runtime object
    const Time& runTime_;

    // Reference to the function object dictionary
    const dictionary& dict_;

public:

    // Construct from time and dictionary
    briscolaFunctionObject
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Disallow default bitwise copy construction
    briscolaFunctionObject(const briscolaFunctionObject&) = delete;

    // Destructor
    virtual ~briscolaFunctionObject();

    // Access to the runtime object
    inline const Time& runTime()
    {
        return runTime_;
    }

    // Access to the function object dictionary
    inline const dictionary& dict()
    {
        return dict_;
    }

    // Return the list of fields required
    virtual wordList fields() const
    {
        return wordList(0);
    }

    // Disallow default bitwise assignment
    void operator=(const briscolaFunctionObject&) = delete;
};

}

}

#endif