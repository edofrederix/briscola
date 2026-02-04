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
                    this->alpha_.surfaceTension()
                )
            );
    }

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
