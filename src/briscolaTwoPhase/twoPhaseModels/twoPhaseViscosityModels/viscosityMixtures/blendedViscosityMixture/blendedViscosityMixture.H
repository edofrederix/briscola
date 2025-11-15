#ifndef blendedViscosityMixture_H
#define blendedViscosityMixture_H

#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class ViscosityModel>
class blendedViscosityMixture
:
    public ViscosityModel
{
protected:

    // Blending factor (defaults to 8)
    const scalar C_;

public:

    // Construct from IOdictionary
    blendedViscosityMixture(const fvMesh&, const IOdictionary&);

    // Copy constructor
    blendedViscosityMixture(const blendedViscosityMixture&);

    // Blending factor
    virtual scalar C() const { return C_; }

    // Correct the two-phase model
    virtual void correct();
};

}

}

}

#ifdef NoRepository
    #include "blendedViscosityMixture.C"
#endif

#endif
