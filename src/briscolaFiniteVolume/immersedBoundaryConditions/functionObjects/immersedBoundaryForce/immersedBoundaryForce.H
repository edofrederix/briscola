#ifndef functionObjects_immersedBoundaryForce_H
#define functionObjects_immersedBoundaryForce_H

#include "briscolaFunctionObject.H"
#include "fvMesh.H"
#include "OFstream.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

template<class MeshType>
class immersedBoundaryForce
:
    public briscolaFunctionObject
{
    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Name of target immersed boundary
    word ibName_;

    // Velocity field name (defaults to U)
    word UName_;

    // Pointer to the immersed boundary
    const immersedBoundary<MeshType>* boundaryPtr_;

    // Calculated force
    List<typename MeshType::vectorType> force_;

    // Output file stream
    autoPtr<OFstream> filePtr_;

public:

    // Runtime type information
    TypeName("immersedBoundaryForce");

    // Construct from time and dictionary
    immersedBoundaryForce
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Disallow default bitwise copy construction
    immersedBoundaryForce(const immersedBoundaryForce&) = delete;

    // Destructor
    virtual ~immersedBoundaryForce();

    // Code executed at the end of each time step loop
    virtual bool execute();

    // Code executed at the end of each time step loop, after execute()
    virtual bool write();

    // Code executed at the end of the simulation
    virtual bool end();

    // Disallow default bitwise assignment
    void operator=(const immersedBoundaryForce&) = delete;
};

}

}

}

}

#endif
