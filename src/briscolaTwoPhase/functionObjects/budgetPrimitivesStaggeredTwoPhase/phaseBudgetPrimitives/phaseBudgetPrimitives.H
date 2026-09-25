#ifndef phaseBudgetPrimitives_H
#define phaseBudgetPrimitives_H

#include "dictionary.H"
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

class budgetPrimitivesStaggeredTwoPhase;

class phaseBudgetPrimitives
{
    const budgetPrimitivesStaggeredTwoPhase& bp_;

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Phase density
    const scalar rho_;

    // Phase viscosity
    const scalar mu_;

    // Phase volume fraction field
    const colocatedScalarField& f_;

    // Colocated velocity field
    const colocatedVectorField& Uc_;

    // Pressure
    const colocatedScalarField& p_;



    // Reconstructed surface tension potential
    const colocatedVectorField& si_;

    // Strain rate tensor
    const colocatedTensorField& Sij_;

    // Velocity gradient
    const colocatedTensorField& djui_;

    // Gradient of phase fraction field
    colocatedVectorField djf_;

    // Phase density ratio
    colocatedScalarField gamma_;

    // Gradient of phase density ratio
    colocatedVectorField djG_;

    // f * ui
    colocatedVectorField fui_;

    // f * ui * ui
    colocatedScalarField fuiui_;

    // djf * ui * ui
    colocatedVectorField djfuiui_;

    // f * ui * uj
    colocatedTensorField fuiuj_;

    // ui * djf
    colocatedTensorField uidjf_;

    // dj (f * ui)
    colocatedTensorField djfui_;

    // f * djui
    colocatedTensorField fdjui_;

    // dj (f * ui * uj)
    colocatedVectorField djfuiuj_;

    // dj (f * ui * ui * uj)
    colocatedScalarField djfuiuiuj_;
    colocatedScalarField djfuiuiujalt_;

    // f * Q
    colocatedScalarField fQ_;

    // dj (f * Q)
    colocatedVectorField djfQ_;

    // f * Q * uj
    colocatedVectorField fQui_;

    // dj (f * Q * uj)
    colocatedTensorField djfQui_;

    // dj (f * nu * Sij) & ui
    colocatedScalarField uidjfnuSij_;

    // dj (f * nu * Sij)
    colocatedScalarField djfnuSxj_;
    colocatedScalarField djfnuSyj_;
    colocatedScalarField djfnuSzj_;

    // f * nu * Sij
    colocatedTensorField fnuSij_;
    colocatedVectorField fnuSxj_;
    colocatedVectorField fnuSyj_;
    colocatedVectorField fnuSzj_;

    // f * nu * Sij * djui
    colocatedScalarField fnuSijdjui_;

    // f * Tij
    colocatedTensorField fTij_;

    // rho * Tij * djG
    colocatedVectorField rhoTijdjG_;

    // rho * Tij * djG * ui
    colocatedScalarField rhoTijuidjG_;

    // f * si
    colocatedVectorField fsi_;

    // f * ui * si
    colocatedScalarField fuisi_;

    // f * f
    colocatedScalarField ff_;

    // Tij
    colocatedTensorField Tij_;

public:

    // Runtime type information
    TypeName("phaseBudgetPrimitives");

    // Construct from time and dictionary
    phaseBudgetPrimitives
    (
        const budgetPrimitivesStaggeredTwoPhase& bpstp,
        const colocatedScalarField& f,
        const label phase
    );

    // Disallow default bitwise copy construction
    phaseBudgetPrimitives(const phaseBudgetPrimitives&) = delete;

    // Destructor
    virtual ~phaseBudgetPrimitives();

    // Compute budget primitives
    void correct();

    // Disallow default bitwise assignment
    void operator=(const phaseBudgetPrimitives&) = delete;
};

}

}

}

}

#endif



