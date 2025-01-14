#ifndef functionObjects_colocatedReconstruction_H
#define functionObjects_colocatedReconstruction_H

#include "briscolaFunctionObject.H"
#include "meshFields.H"
#include "exSchemes.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

class colocatedReconstruction
:
    public briscolaFunctionObject
{
    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Fields to be reconstructed
    wordList fields_;

    // Reconstructed vectors
    PtrList<colocatedVectorField> vectorFields_;

    // Reconstructed tensors
    PtrList<colocatedTensorField> tensorFields_;

    // Initialize
    void init();

public:

    // Runtime type information
    TypeName("colocatedReconstruction");

    // Construct from time and dictionary
    colocatedReconstruction
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Disallow default bitwise copy construction
    colocatedReconstruction(const colocatedReconstruction&) = delete;

    // Destructor
    virtual ~colocatedReconstruction();

    // Code executed at the end of each time step loop
    virtual bool execute();

    // Code executed at the end of each time step loop, after execute()
    virtual bool write();

    // Code executed at the end of the simulation
    virtual bool end();

    // Disallow default bitwise assignment
    void operator=(const colocatedReconstruction&) = delete;
};

}

}

}

}

#endif
