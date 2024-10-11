#ifndef coloBlendedViscosityRestrictionScheme_H
#define coloBlendedViscosityRestrictionScheme_H

#include "blendedViscosityRestrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class coloBlendedViscosityRestrictionScheme
:
    public blendedViscosityRestrictionScheme<colocated>
{
public:

    // Runtime type information
    TypeName("blendedViscosity");

    // Construct from stream and mesh
    coloBlendedViscosityRestrictionScheme
    (
        const fvMesh&,
        Istream& = scheme::nullStream
    );

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

#endif
