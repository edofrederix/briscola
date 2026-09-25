#ifndef functionObjects_budgetPrimitivesStaggeredTwoPhase_H
#define functionObjects_budgetPrimitivesStaggeredTwoPhase_H

#include "briscolaFunctionObject.H"
#include "IOdictionary.H"
#include "meshFields.H"
#include "exSchemes.H"
#include "phaseBudgetPrimitives.H"
#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

class budgetPrimitivesStaggeredTwoPhase
:
    public briscolaFunctionObject
{
    friend class phaseBudgetPrimitives;

    // Number of timesteps per budget calculation iteration
    int samplingSteps_;

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Two-phase model
    const twoPhaseModel& tpm_;

    // Colocated velocity field
    const colocatedVectorField& Uc_;

    // Pressure
    const colocatedScalarField& p_;

    // 1 - alpha
    colocatedScalarField f1_;

    // alpha
    const colocatedScalarField& f2_;

    // Density of phase 1
    const scalar rho1_;

    // Density of phase 2
    const scalar rho2_;

    // Viscosity of phase 1
    const scalar mu1_;

    // Viscosity of phase 2
    const scalar mu2_;

    // Mixture density
    colocatedScalarField rho_;

    // Mixture viscosity
    colocatedScalarField mu_;

    // Reconstructed surface tension potential
    colocatedVectorField si_;

    // Gradient of velocity
    colocatedTensorField djui_;

    // Strain rate tensor
    colocatedTensorField Sij_;

    // Budget primitives for phase 1
    phaseBudgetPrimitives bp1_;

    // Budget primitives for phase 2
    phaseBudgetPrimitives bp2_;


public:

    // Runtime type information
    TypeName("budgetPrimitivesStaggeredTwoPhase");

    // Construct from time and dictionary
    budgetPrimitivesStaggeredTwoPhase
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Destructor
    virtual ~budgetPrimitivesStaggeredTwoPhase();

    // Code executed at the end of each time step loop
    virtual bool execute();

    // Code executed at the end of each time step loop, after execute()
    virtual bool write();

    // Code executed at the end of the simulation
    virtual bool end();

    // Disallow default bitwise assignment
    void operator=(const budgetPrimitivesStaggeredTwoPhase&) = delete;
};

}

}

}

}

#endif



