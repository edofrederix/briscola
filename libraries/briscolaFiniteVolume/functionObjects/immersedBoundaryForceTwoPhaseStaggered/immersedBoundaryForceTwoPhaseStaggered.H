#ifndef functionObjects_immersedBoundaryForceTwoPhaseStaggered_H
#define functionObjects_immersedBoundaryForceTwoPhaseStaggered_H

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

class immersedBoundaryForceTwoPhaseStaggered
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

    // Number of target immersed boundary
    label number_;

public:

    // Runtime type information
    TypeName("immersedBoundaryForceTwoPhaseStaggered");

    // Construct from time and dictionary
    immersedBoundaryForceTwoPhaseStaggered
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Disallow default bitwise copy construction
    immersedBoundaryForceTwoPhaseStaggered(const immersedBoundaryForceTwoPhaseStaggered&) = delete;

    // Destructor
    virtual ~immersedBoundaryForceTwoPhaseStaggered();

    // Code executed at the end of each time step loop
    virtual bool execute();

    // Code executed at the end of each time step loop, after execute()
    virtual bool write();

    // Code executed at the end of the simulation
    virtual bool end();

    // Disallow default bitwise assignment
    void operator=(const immersedBoundaryForceTwoPhaseStaggered&) = delete;
};

}

}

}

}

#endif
