#ifndef coloFluxRestrictionScheme_H
#define coloFluxRestrictionScheme_H

#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class coloFluxRestrictionScheme
:
    public restrictionScheme<lowerFaceScalar,colocated>
{
public:

    // Runtime type information
    TypeName("flux");

    // Construct from dictionary and mesh
    coloFluxRestrictionScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    coloFluxRestrictionScheme(const fvMesh&);

    // In-place restriction from a fine cell level to an adjacent coarse cell
    // level
    virtual void restrict
    (
        meshDirection<lowerFaceScalar,colocated>&,
        const meshDirection<lowerFaceScalar,colocated>&,
        const bool = false
    );
};

}

}

}

#ifdef NoRepository
    #include "coloFluxRestrictionScheme.C"
#endif

#endif
