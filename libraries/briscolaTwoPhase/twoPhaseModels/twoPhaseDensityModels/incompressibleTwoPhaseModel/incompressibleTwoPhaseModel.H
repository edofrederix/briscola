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
    const scalar Rho1_;
    const scalar Rho2_;

private:

    // Disallow default bitwise assignment
    void operator=(const incompressibleTwoPhaseModel&);

public:

    // Runtime type information
    TypeName("incompressibleTwoPhaseModel");

    // Construct from IOdictionary
    incompressibleTwoPhaseModel(const fvMesh&, const IOdictionary&);

    // Copy constructor
    incompressibleTwoPhaseModel(const incompressibleTwoPhaseModel&);

    // Destructor
    virtual ~incompressibleTwoPhaseModel();

    // Mass density of the first phase
    virtual scalar Rho1() { return Rho1_; }

    // Mass density of the second phase
    virtual scalar Rho2() { return Rho2_; }

    // Correct the two-phase mixture
    virtual void correctMixture();
};

}

}

}

#endif
