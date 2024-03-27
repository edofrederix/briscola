#ifndef functionObjects_budgetPrimitivesStaggered_H
#define functionObjects_budgetPrimitivesStaggered_H

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

class budgetPrimitivesStaggered
:
    public briscolaFunctionObject
{
    // Reference to the mesh
    const fvMesh& fvMsh_;

    // velocity squared
    staggeredScalarField U2_;

    // p*p
    colocatedScalarField p2_;

    // Reynolds stress tensor
    // Based on reconstructed colocated velocity,
    // so turn on colocatedReconstruction switch.
    colocatedSymmTensorField Rij_;

    // u_i*u_i*u_k
    colocatedTensorField uiuiuk_;

    // p*ui
    colocatedVectorField pui_;

    // Gradient of velocity
    colocatedTensorField djui_;

    // Strain rate tensor
    colocatedTensorField Sij_;

    // Strain rate tensor squared
    colocatedTensorField eij_;

public:

    // Runtime type information
    TypeName("budgetPrimitivesStaggered");

    // Construct from time and dictionary
    budgetPrimitivesStaggered
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Disallow default bitwise copy construction
    budgetPrimitivesStaggered(const budgetPrimitivesStaggered&) = delete;

    // Destructor
    virtual ~budgetPrimitivesStaggered();

    // Code executed at the end of each time step loop
    virtual bool execute();

    // Code executed at the end of each time step loop, after execute()
    virtual bool write();

    // Code executed at the end of the simulation
    virtual bool end();

    // Disallow default bitwise assignment
    void operator=(const budgetPrimitivesStaggered&) = delete;
};

}

}

}

}

#endif
