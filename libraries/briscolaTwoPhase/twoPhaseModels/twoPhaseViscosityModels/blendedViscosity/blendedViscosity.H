#ifndef blendedViscosity_H
#define blendedViscosity_H

#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
class blendedViscosity
:
    public BaseModel
{
protected:

    // The dynamic viscosities of the two phases
    const scalar mu1_;
    const scalar mu2_;

    // Blending factor (defaults to 4)
    const scalar C_;

private:

    // Disallow default bitwise assignment
    void operator=(const blendedViscosity&);

public:

    // Return the instantiated type name
    static word typeName()
    {
        return "blendedViscosity<" + word(BaseModel::typeName) + '>';
    }

    // Construct from IOdictionary
    blendedViscosity(const fvMesh&, const IOdictionary&);

    // Copy constructor
    blendedViscosity(const blendedViscosity&);

    // Destructor
    virtual ~blendedViscosity();

    // Dynamic viscosity of the first phase
    virtual scalar mu1() const { return mu1_; }

    // Dynamic viscosity of the second phase
    virtual scalar mu2() const { return mu2_; }

    // Blending factor
    virtual scalar C() const { return C_; }

    // Correct the two-phase mixture
    virtual void correctMixture();
};

}

}

}

#ifdef NoRepository
    #include "blendedViscosity.C"
#endif

#endif
