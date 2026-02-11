#ifndef functionObjects_bubblePipe2D_H
#define functionObjects_bubblePipe2D_H

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

class bubblePipe2D
:
    public briscolaFunctionObject
{
    // Reference to the mesh
    const fvMesh& fvMsh_;

public:

    // Runtime type information
    TypeName("bubblePipe2D");

    // Construct from time and dictionary
    bubblePipe2D
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Disallow default bitwise copy construction
    bubblePipe2D(const bubblePipe2D&) = delete;

    // Destructor
    virtual ~bubblePipe2D();

    // Read the dictionary. Called at construction.
    virtual bool read(const dictionary&);

    // Code executed at the end of each time step loop
    virtual bool execute() { return true; }

    // Code executed at the end of each time step loop, after execute()
    virtual bool write() { return true; }

    // Code executed at the end of the simulation
    virtual bool end() { return true; }

    // Disallow default bitwise assignment
    void operator=(const bubblePipe2D&) = delete;
};

}

}

}

}

#endif
