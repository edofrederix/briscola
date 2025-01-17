#ifndef blendedViscosityMixture_H
#define blendedViscosityMixture_H

#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
class blendedViscosityMixture
:
    public BaseModel
{
protected:

    // Blending factor (defaults to 8)
    const scalar C_;

private:

    // Disallow default bitwise assignment
    void operator=(const blendedViscosityMixture&);

public:

    // Return the instantiated type name
    static word typeName()
    {
        return "blendedViscosityMixture<" + word(BaseModel::typeName) + '>';
    }

    // Construct from IOdictionary
    blendedViscosityMixture(const fvMesh&, const IOdictionary&);

    // Copy constructor
    blendedViscosityMixture(const blendedViscosityMixture&);

    // Destructor
    virtual ~blendedViscosityMixture();

    // Blending factor
    virtual scalar C() const { return C_; }

    // Correct the two-phase mixture
    virtual void correctMixture();
};

}

}

}

#ifdef NoRepository
    #include "blendedViscosityMixture.C"
#endif

#endif
