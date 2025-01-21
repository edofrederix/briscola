#ifndef incompressibleTwoPhaseModel_H
#define incompressibleTwoPhaseModel_H

#include "TwoPhaseModel.H"
#include "vof.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
class incompressibleTwoPhaseModel
:
    public TwoPhaseModel<MeshType>
{
protected:

    // The mass densities of the two phases
    const scalar rho1_;
    const scalar rho2_;

private:

    // Disallow default bitwise assignment
    void operator=(const incompressibleTwoPhaseModel&);

public:

    // Construct from IOdictionary
    incompressibleTwoPhaseModel(const fvMesh&, const IOdictionary&);

    // Copy constructor
    incompressibleTwoPhaseModel(const incompressibleTwoPhaseModel&);

    // Destructor
    virtual ~incompressibleTwoPhaseModel();

    // Mass density of the first phase
    virtual scalar rho1() { return rho1_; }

    // Mass density of the second phase
    virtual scalar rho2() { return rho2_; }

    // Return the colocated mass density
    virtual tmp<colocatedScalarField> coloRho() const;

    // Correct the two-phase mixture
    virtual void correctMixture();
};

}

}

}

#ifdef NoRepository
    #include "incompressibleTwoPhaseModel.C"
#endif

#endif
