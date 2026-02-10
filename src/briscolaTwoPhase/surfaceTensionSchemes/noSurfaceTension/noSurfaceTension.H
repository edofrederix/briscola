#ifndef noSurfaceTension_H
#define noSurfaceTension_H

#include "surfaceTensionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SigmaModel>
class noSurfaceTension
:
    public SigmaModel
{
public:

    // Runtime type information
    TypeName("none");

    // Construct from components
    noSurfaceTension
    (
        const fvMesh&,
        const dictionary&,
        const normalScheme&,
        const colocatedScalarField&
    );

    // Copy constructor
    noSurfaceTension(const noSurfaceTension&);

    // Destructor
    virtual ~noSurfaceTension();

    // Clone surfaceTensionScheme
    virtual autoPtr<surfaceTensionScheme> clone() const
    {
        return autoPtr<surfaceTensionScheme>(new noSurfaceTension(*this));
    }

    // Correct the surface tension force
    virtual void correct();
};

}

}

}

#ifdef NoRepository
    #include "noSurfaceTension.C"
#endif

#endif
