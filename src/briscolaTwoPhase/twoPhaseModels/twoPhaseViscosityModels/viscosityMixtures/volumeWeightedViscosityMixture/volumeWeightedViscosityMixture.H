#ifndef volumeWeightedViscosityMixture_H
#define volumeWeightedViscosityMixture_H

#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
class volumeWeightedViscosityMixture
:
    public BaseModel
{
private:

    // Disallow default bitwise assignment
    void operator=(const volumeWeightedViscosityMixture&);

public:

    // Construct from IOdictionary
    volumeWeightedViscosityMixture(const fvMesh&, const IOdictionary&);

    // Copy constructor
    volumeWeightedViscosityMixture(const volumeWeightedViscosityMixture&);

    // Destructor
    virtual ~volumeWeightedViscosityMixture();

    // Correct the two-phase mixture
    virtual void correctMixture();
};

}

}

}

#ifdef NoRepository
    #include "volumeWeightedViscosityMixture.C"
#endif

#endif
