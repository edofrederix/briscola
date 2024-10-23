#ifndef stagBlendedViscosityMixtureRestrictionScheme_H
#define stagBlendedViscosityMixtureRestrictionScheme_H

#include "blendedViscosityMixtureRestrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class stagBlendedViscosityMixtureRestrictionScheme
:
    public blendedViscosityMixtureRestrictionScheme<staggered>
{
public:

    // Runtime type information
    TypeName("blendedViscosityMixture");

    // Construct from stream and mesh
    stagBlendedViscosityMixtureRestrictionScheme
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
