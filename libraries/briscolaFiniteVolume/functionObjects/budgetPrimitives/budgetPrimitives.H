#ifndef functionObjects_budgetPrimitives_H
#define functionObjects_budgetPrimitives_H

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

class budgetPrimitives
:
    public briscolaFunctionObject
{
    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Name of colocated velocity field
    word uName_;

    // u*u
    colocatedScalarField u2_;

    // v*v
    colocatedScalarField v2_;

    // w*w
    colocatedScalarField w2_;

public:

    // Runtime type information
    TypeName("budgetPrimitives");

    // Construct from time and dictionary
    budgetPrimitives
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Disallow default bitwise copy construction
    budgetPrimitives(const budgetPrimitives&) = delete;

    // Destructor
    virtual ~budgetPrimitives();

    // Code executed at the end of each time step loop
    virtual bool execute();

    // Code executed at the end of each time step loop, after execute()
    virtual bool write();

    // Code executed at the end of the simulation
    virtual bool end();

    // Disallow default bitwise assignment
    void operator=(const budgetPrimitives&) = delete;
};

}

}

}

}

#endif
