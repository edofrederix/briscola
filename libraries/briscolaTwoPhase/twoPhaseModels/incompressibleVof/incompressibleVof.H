#ifndef incompressibleVof_H
#define incompressibleVof_H

#include "twoPhaseModel.H"
#include "vof.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class incompressibleVof
:
    public twoPhaseModel
{
protected:

    // The Vof model
    autoPtr<vof> vfPtr_;

    // The mass densities of the two components
    const scalar rho1_;
    const scalar rho2_;

    // The dynamic viscosities of the two components
    const scalar mu1_;
    const scalar mu2_;

    // Name of the colocated volumetric mixture flux field (defaults to phi)
    const word fluxName_;

    // Name of the staggered velocity field, only needed when the colocated
    // volumetric mixture flux field is not found (defaults to U)
    const word velocityName_;

private:

    // Disallow default bitwise assignment
    void operator=(const incompressibleVof&);

public:

    // Runtime type information
    TypeName("incompressibleVof");

    // Construct from IOdictionary
    incompressibleVof(const IOdictionary&, const fvMesh&);

    // Copy constructor
    incompressibleVof(const incompressibleVof&);

    // Clone incompressibleVof
    virtual autoPtr<twoPhaseModel> clone() const
    {
        return autoPtr<twoPhaseModel>(new incompressibleVof(*this));
    }

    // Destructor
    ~incompressibleVof();

    // Return access to the Vof model
    inline vof& vf() { return vfPtr_(); }

    // Return const access to the Vof model
    inline const vof& vf() const { return vfPtr_; }

    // Correct the two-phase model
    void correct();
};

}

}

}

#endif
