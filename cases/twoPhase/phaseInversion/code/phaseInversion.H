#ifndef functionObjects_phaseInversion_H
#define functionObjects_phaseInversion_H

#include "briscolaFunctionObject.H"
#include "OFstream.H"
#include "IOdictionary.H"
#include "meshFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

class phaseInversion
:
    public briscolaFunctionObject
{
    // Output file stream
    autoPtr<OFstream> filePtr_;

    // Densities
    scalar rho1_;
    scalar rho2_;

    // Gravity
    scalar g_;

public:

    // Runtime type information
    TypeName("phaseInversion");

    // Construct from time and dictionary
    phaseInversion
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Read the dictionary. Called at construction.
    virtual bool read(const dictionary&);

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