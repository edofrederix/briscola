#ifndef NewtonianViscosity_H
#define NewtonianViscosity_H

#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
class NewtonianViscosity
:
    public BaseModel
{
protected:

    // The dynamic viscosities of the two phases
    const scalar mu1_;
    const scalar mu2_;

private:

    // Disallow default bitwise assignment
    void operator=(const NewtonianViscosity&);

public:

    // Return the instantiated type name
    static word typeName()
    {
        return "NewtonianViscosity<" + word(BaseModel::typeName) + '>';
    }

    // Construct from IOdictionary
    NewtonianViscosity(const fvMesh&, const IOdictionary&);

    // Copy constructor
    NewtonianViscosity(const NewtonianViscosity&);

    // Destructor
    virtual ~NewtonianViscosity();

    // Dynamic viscosity of the first phase
    virtual scalar mu1() const { return mu1_; }

    // Dynamic viscosity of the second phase
    virtual scalar mu2() const { return mu2_; }

    // Correct the two-phase mixture
    virtual void correctMixture();
};

}

}

}

#ifdef NoRepository
    #include "NewtonianViscosity.C"
#endif

#endif
