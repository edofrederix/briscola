#ifndef stagBlendedViscosityRestrictionScheme_H
#define stagBlendedViscosityRestrictionScheme_H

#include "blendedViscosityRestrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class stagBlendedViscosityRestrictionScheme
:
    public blendedViscosityRestrictionScheme<staggered>
{
public:

    // Runtime type information
    TypeName("blendedViscosity");

    // Construct from stream and mesh
    stagBlendedViscosityRestrictionScheme
    (
        const fvMesh&,
        Istream& = scheme::nullStream
    );

    // In-place restriction from a fine cell level to an adjacent coarse cell
    // level
    virtual void restrict
    (
        meshDirection<lowerFaceScalar,staggered>&,
        const meshDirection<lowerFaceScalar,staggered>&,
        const bool = false
    );
};

}

}

}

#endif
