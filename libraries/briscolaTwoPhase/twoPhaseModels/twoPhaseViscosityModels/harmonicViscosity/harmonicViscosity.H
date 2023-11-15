#ifndef harmonicViscosity_H
#define harmonicViscosity_H

#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
class harmonicViscosity
:
    public BaseModel
{
protected:

    // The dynamic viscosities of the two phases
    const scalar mu1_;
    const scalar mu2_;

private:

    // Disallow default bitwise assignment
    void operator=(const harmonicViscosity&);

public:

    // Return the instantiated type name
    static word typeName()
    {
        return "harmonicViscosity<" + word(BaseModel::typeName) + '>';
    }

    // Construct from IOdictionary
    harmonicViscosity(const fvMesh&, const IOdictionary&);

    // Copy constructor
    harmonicViscosity(const harmonicViscosity&);

    // Destructor
    virtual ~harmonicViscosity();

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
    #include "harmonicViscosity.C"
#endif

#endif
