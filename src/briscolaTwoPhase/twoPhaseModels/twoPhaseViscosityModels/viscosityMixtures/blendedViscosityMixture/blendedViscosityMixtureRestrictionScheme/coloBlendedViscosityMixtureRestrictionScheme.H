#ifndef coloBlendedViscosityMixtureRestrictionScheme_H
#define coloBlendedViscosityMixtureRestrictionScheme_H

#include "blendedViscosityMixtureRestrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class coloBlendedViscosityMixtureRestrictionScheme
:
    public blendedViscosityMixtureRestrictionScheme<colocated>
{
public:

    // Runtime type information
    TypeName("blendedViscosityMixture");

    // Construct from stream and mesh
    coloBlendedViscosityMixtureRestrictionScheme
    (
        const fvMesh&,
        Istream& = scheme::nullStream
    );

    // In-place restriction from a fine cell level to an adjacent coarse cell
    // level
    virtual void restrict
    (
        meshDirection<faceScalar,colocated>&,
        const meshDirection<faceScalar,colocated>&,
        const bool = false
    );
};

}

}

}

#endif
