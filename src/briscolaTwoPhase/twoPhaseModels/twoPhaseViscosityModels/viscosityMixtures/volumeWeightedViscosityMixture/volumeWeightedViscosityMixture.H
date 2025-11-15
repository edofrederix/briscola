#ifndef volumeWeightedViscosityMixture_H
#define volumeWeightedViscosityMixture_H

#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class ViscosityModel>
class volumeWeightedViscosityMixture
:
    public ViscosityModel
{
public:

    // Construct from IOdictionary
    volumeWeightedViscosityMixture(const fvMesh&, const IOdictionary&);

    // Copy constructor
    volumeWeightedViscosityMixture(const volumeWeightedViscosityMixture&);

    // Correct the two-phase model
    virtual void correct();
};

}

}

}

#ifdef NoRepository
    #include "volumeWeightedViscosityMixture.C"
#endif

#endif
