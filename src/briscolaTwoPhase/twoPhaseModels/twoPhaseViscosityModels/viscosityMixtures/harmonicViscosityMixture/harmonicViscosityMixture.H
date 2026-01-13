#ifndef harmonicViscosityMixture_H
#define harmonicViscosityMixture_H

#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class ViscosityModel>
class harmonicViscosityMixture
:
    public ViscosityModel
{
public:

    // Construct from IOdictionary
    harmonicViscosityMixture(const fvMesh&, const IOdictionary&);

    // Copy constructor
    harmonicViscosityMixture(const harmonicViscosityMixture&);

    // Correct the two-phase model
    virtual void correct();
};

}

}

}

#ifdef NoRepository
    #include "harmonicViscosityMixture.C"
#endif

#endif
