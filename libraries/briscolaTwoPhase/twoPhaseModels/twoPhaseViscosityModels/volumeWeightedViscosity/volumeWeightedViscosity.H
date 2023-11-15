#ifndef volumeWeightedViscosity_H
#define volumeWeightedViscosity_H

#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
class volumeWeightedViscosity
:
    public BaseModel
{
protected:

    // The dynamic viscosities of the two phases
    const scalar mu1_;
    const scalar mu2_;

private:

    // Disallow default bitwise assignment
    void operator=(const volumeWeightedViscosity&);

public:

    // Construct from IOdictionary
    volumeWeightedViscosity(const fvMesh&, const IOdictionary&);

    // Copy constructor
    volumeWeightedViscosity(const volumeWeightedViscosity&);

    // Destructor
    virtual ~volumeWeightedViscosity();

    // Dynamic viscosity of the first phase
    virtual scalar mu1() { return mu1_; }

    // Dynamic viscosity of the second phase
    virtual scalar mu2() { return mu2_; }

    // Correct the two-phase mixture
    virtual void correctMixture();
};

}

}

}

#ifdef NoRepository
    #include "volumeWeightedViscosity.C"
#endif

#endif
