#ifndef Brackbill_H
#define Brackbill_H

#include "surfaceTensionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SigmaModel>
class Brackbill
:
    public SigmaModel
{
public:

    // Runtime type information
    TypeName("Brackbill");

    // Construct from components
    Brackbill
    (
        const twoPhaseModel&,
        const dictionary&
    );

    // Copy constructor
    Brackbill(const Brackbill&);

    // Destructor
    virtual ~Brackbill();

    // Clone surfaceTensionScheme
    virtual autoPtr<surfaceTensionScheme> clone() const
    {
        return autoPtr<surfaceTensionScheme>(new Brackbill(*this));
    }

    // Correct the surface tension force
    virtual void correct();
};

}

}

}

#ifdef NoRepository
    #include "Brackbill.C"
#endif

#endif
