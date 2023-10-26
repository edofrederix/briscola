#ifndef geometricVof_H
#define geometricVof_H

#include "vof.H"
#include "LVE.H"
#include "normalScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class geometricVof
:
    public vof
{
protected:

    // Rectilinear bool (stored here for efficiency)
    bool rectilinear_;

    // Local volume enforcement problem solver
    LVE lve_;

    // Normal scheme
    autoPtr<normalScheme> normalSchemePtr_;

private:

    // Disallow default bitwise assignment
    void operator=(const geometricVof&);

public:

    // Runtime type information
    TypeName("geometricVof");

    // Construct from dictionary
    geometricVof(const dictionary&, const fvMesh&);

    // Copy constructor
    geometricVof(const geometricVof&);

    // Destructor
    ~geometricVof();

    // Reference to the LVE solver
    inline const LVE& lve() const { return lve_; }

    // Reference to the normal scheme
    inline normalScheme& normal() { return normalSchemePtr_(); }
};

}

}

}

#endif
