#ifndef twoPhaseVof_H
#define twoPhaseVof_H

#include "vof.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class ViscosityModel>
class twoPhaseVof
:
    public ViscosityModel
{
protected:

    // Normal scheme
    autoPtr<normalScheme> normalSchemePtr_;

    // Surface tension scheme
    autoPtr<surfaceTensionScheme> surfaceTensionSchemePtr_;

    // The Vof model
    autoPtr<vof> vfPtr_;

public:

    // Runtime type information
    TypeName("twoPhaseVof");

    // Construct from IOdictionary
    twoPhaseVof(const fvMesh&, const IOdictionary&);

    // Copy constructor
    twoPhaseVof(const twoPhaseVof&);

    // Clone twoPhaseVof
    virtual autoPtr<twoPhaseModel> clone() const
    {
        return autoPtr<twoPhaseModel>(new twoPhaseVof(*this));
    }

    // Return surface tension flux
    virtual tmp<colocatedScalarFaceField> flux()
    {
        return tmp<colocatedScalarFaceField>
            (
                static_cast<colocatedScalarFaceField&>
                (
                    surfaceTensionSchemePtr_()
                )
            );
    }

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
    #include "twoPhaseVof.C"
#endif

#endif
