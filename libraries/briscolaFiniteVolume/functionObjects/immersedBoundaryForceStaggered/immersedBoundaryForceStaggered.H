#ifndef functionObjects_immersedBoundaryForceStaggered_H
#define functionObjects_immersedBoundaryForceStaggered_H

#include "briscolaFunctionObject.H"
#include "fvMesh.H"
#include "vector.H"
#include "meshField.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

class immersedBoundaryForceStaggered
:
    public briscolaFunctionObject
{

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Volumetric force
    staggeredScalarField volForce_;

    // Name of target immersed boundary
    word name_;

    // Immersed boundary method
    word method_;

    // Force from explicit source (assumed constant in space and time)
    vector exSourceForce_;

    // Number of target immersed boundary
    label number_;

public:

    // Runtime type information
    TypeName("immersedBoundaryForceStaggered");

    // Construct from time and dictionary
    immersedBoundaryForceStaggered
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Disallow default bitwise copy construction
    immersedBoundaryForceStaggered(const immersedBoundaryForceStaggered&) = delete;

    // Destructor
    virtual ~immersedBoundaryForceStaggered();

    // Code executed at the end of each time step loop
    virtual bool execute();

    // Code executed at the end of each time step loop, after execute()
    virtual bool write();

    // Code executed at the end of the simulation
    virtual bool end();

    // Compute force from explicit source term (assumed constant in space and time)
    void computeExplicitSourceForce();

    // Disallow default bitwise assignment
    void operator=(const immersedBoundaryForceStaggered&) = delete;
};

}

}

}

}

#endif
