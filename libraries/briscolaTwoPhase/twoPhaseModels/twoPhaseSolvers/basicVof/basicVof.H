#ifndef basicVof_H
#define basicVof_H

#include "vof.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
class basicVof
:
    public BaseModel
{
protected:

    // The Vof model
    autoPtr<vof> vfPtr_;

    // Name of the colocated volumetric mixture flux field (defaults to phi)
    const word fluxName_;

    // Name of the staggered velocity field, only needed when the colocated
    // volumetric mixture flux field is not found (defaults to U)
    const word velocityName_;

private:

    // Disallow default bitwise assignment
    void operator=(const basicVof&);

public:

    // Runtime type information
    TypeName("basicVof");

    // Construct from IOdictionary
    basicVof(const IOdictionary&, const fvMesh&);

    // Copy constructor
    basicVof(const basicVof&);

    // Clone basicVof
    virtual autoPtr<twoPhaseModel> clone() const
    {
        return autoPtr<twoPhaseModel>(new basicVof(*this));
    }

    // Destructor
    virtual ~basicVof();

    // Return access to the Vof model
    inline vof& vf() { return vfPtr_(); }

    // Return const access to the Vof model
    inline const vof& vf() const { return vfPtr_; }

    // Correct the two-phase model
    virtual void correct();
};

}

}

}

#ifdef NoRepository
    #include "basicVof.C"
#endif

#endif
