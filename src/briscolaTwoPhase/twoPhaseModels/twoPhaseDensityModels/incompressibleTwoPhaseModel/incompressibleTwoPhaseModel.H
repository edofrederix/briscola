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

    // Correct the mixture mass density
    virtual void correctRho();

    // Correct the mean mixture mass density
    virtual void correctMeanRho();

public:

    // Construct from IOdictionary
    incompressibleTwoPhaseModel(const fvMesh&, const IOdictionary&);

    // Copy constructor
    incompressibleTwoPhaseModel(const incompressibleTwoPhaseModel&);

    // Mass density of the first phase
    virtual scalar rho1() const { return rho1_; }

    // Mass density of the second phase
    virtual scalar rho2() const { return rho2_; }

    // Correct the two-phase model
    virtual void correct();
};

}

}

}

#ifdef NoRepository
    #include "incompressibleTwoPhaseModel.C"
#endif

#endif
