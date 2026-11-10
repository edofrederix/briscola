#ifndef functionObjects_PILULE_H
#define functionObjects_PILULE_H

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

class PILULE
:
    public briscolaFunctionObject
{
    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Output file stream
    autoPtr<OFstream> filePtr_;

    // Total fluid volume
    scalar fluidVol_;

    // Volume lost to IBM
    scalar removedVol_;

    // Volume lost field
    colocatedScalarField volLossIBM_;

    // Phase 1 velocity
    colocatedVectorField U1_;

    // Phase 1 velocity squared
    colocatedVectorField Usqr1_;

    // Phase 2 velocity
    colocatedVectorField U2_;

    // Phase 2 velocity squared
    colocatedVectorField Usqr2_;

    // Dictionary to communicate missing volume
    IOdictionary piluleDict_;

    // PI controller quantities

        // Error in previous time step
        scalar oldError_;

        // Proportional gain
        scalar kP_;

        // Integral gain
        scalar kI_;

        // Target J2
        scalar targetJ2_;

        // Current forcing value
        scalar forcing_;

public:

    // Runtime type information
    TypeName("PILULE");

    // Construct from time and dictionary
    PILULE
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Disallow default bitwise copy construction
    PILULE(const PILULE&) = delete;

    // Destructor
    virtual ~PILULE();

    // Read the dictionary. Called at construction.
    virtual bool read(const dictionary&);

    // Code executed at the end of each time step loop
    virtual bool execute();

    // Code executed at the end of each time step loop, after execute()
    virtual bool write() { return true; }

    // Code executed at the end of the simulation
    virtual bool end() { return true; }

    // Disallow default bitwise assignment
    void operator=(const PILULE&) = delete;
};

}

}

}

}

#endif
