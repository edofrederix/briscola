#ifndef harmonicViscosityMixture_H
#define harmonicViscosityMixture_H

#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
class harmonicViscosityMixture
:
    public BaseModel
{
private:

    // Disallow default bitwise assignment
    void operator=(const harmonicViscosityMixture&);

public:

    // Return the instantiated type name
    static word typeName()
    {
        return "harmonicViscosityMixture<" + word(BaseModel::typeName) + '>';
    }

    // Construct from IOdictionary
    harmonicViscosityMixture(const fvMesh&, const IOdictionary&);

    // Copy constructor
    harmonicViscosityMixture(const harmonicViscosityMixture&);

    // Destructor
    virtual ~harmonicViscosityMixture();

    // Correct the two-phase mixture
    virtual void correctMixture();
};

}

}

}

#ifdef NoRepository
    #include "harmonicViscosityMixture.C"
#endif

#endif
