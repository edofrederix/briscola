#ifndef blendedViscosityMixtureRestrictionScheme_H
#define blendedViscosityMixtureRestrictionScheme_H

#include "blendedViscosityMixtureRestrictionSchemeBase.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Forward declaration

template<class MeshType>
class blendedViscosityMixtureRestrictionScheme;


// Colocated

template<>
class blendedViscosityMixtureRestrictionScheme<colocated>
:
    public blendedViscosityMixtureRestrictionSchemeBase<colocated>
{
public:

    // Runtime type information
    TypeName("blendedViscosityMixture");

    // Use base constructors
    using blendedViscosityMixtureRestrictionSchemeBase<colocated>::
        blendedViscosityMixtureRestrictionSchemeBase;

    // // Destructor
    // virtual ~blendedViscosityMixtureRestrictionScheme();

    // In-place restriction from a fine cell level to an adjacent coarse cell
    // level
    virtual void restrict
    (
        meshDirection<faceScalar,colocated>&,
        const meshDirection<faceScalar,colocated>&,
        const bool = false
    );
};


// Staggered

template<>
class blendedViscosityMixtureRestrictionScheme<staggered>
:
    public blendedViscosityMixtureRestrictionSchemeBase<staggered>
{
public:

    // Runtime type information
    TypeName("blendedViscosityMixture");

    // Use base constructors
    using blendedViscosityMixtureRestrictionSchemeBase<staggered>::
        blendedViscosityMixtureRestrictionSchemeBase;

    // // Destructor
    // virtual ~blendedViscosityMixtureRestrictionScheme();

    // In-place restriction from a fine cell level to an adjacent coarse cell
    // level
    virtual void restrict
    (
        meshDirection<faceScalar,staggered>&,
        const meshDirection<faceScalar,staggered>&,
        const bool = false
    );
};

}

}

}

#ifdef NoRepository
    #include "blendedViscosityMixtureRestrictionScheme.C"
#endif

#endif
