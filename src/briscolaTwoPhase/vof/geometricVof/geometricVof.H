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

public:

    // Construct from components
    geometricVof
    (
        const fvMesh&,
        const dictionary&,
        normalScheme&,
        vofField&
    );

    // Copy constructor
    geometricVof(const geometricVof&);

    // Reference to the LVE solver
    inline const LVE& lve() const { return lve_; }
};

}

}

}

#endif
