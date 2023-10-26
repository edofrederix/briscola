#ifndef incompressibleTwoPhaseModel_H
#define incompressibleTwoPhaseModel_H

#include "twoPhaseModel.H"
#include "vof.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class incompressibleTwoPhaseModel
:
    public twoPhaseModel
{
protected:

    // The mass densities of the two phases
    const scalar rho1_;
    const scalar rho2_;

private:

    // Disallow default bitwise assignment
    void operator=(const incompressibleTwoPhaseModel&);

public:

    // Runtime type information
    TypeName("incompressibleTwoPhaseModel");

    // Construct from IOdictionary
    incompressibleTwoPhaseModel(const IOdictionary&, const fvMesh&);

    // Copy constructor
    incompressibleTwoPhaseModel(const incompressibleTwoPhaseModel&);

    // Destructor
    virtual ~incompressibleTwoPhaseModel();

    // Mass density of the first phase
    virtual scalar rho1() { return rho1_; }

    // Mass density of the second phase
    virtual scalar rho2() { return rho2_; }

    // Correct the two-phase mixture
    virtual void correctMixture(const colocatedScalarField&);
};

}

}

}

#endif
