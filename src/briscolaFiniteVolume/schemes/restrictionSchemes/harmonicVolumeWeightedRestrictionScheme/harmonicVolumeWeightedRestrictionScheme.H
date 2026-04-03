#ifndef harmonicVolumeWeightedRestrictionScheme_H
#define harmonicVolumeWeightedRestrictionScheme_H

#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class harmonicVolumeWeightedRestrictionScheme
:
    public restrictionScheme<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("harmonicVolumeWeighted");

    // Use base constructors
    using restrictionScheme<Type,MeshType>::restrictionScheme;

    // In-place restriction from a fine cell level to an adjacent coarse cell
    // level, with option to scale by cell volume
    virtual void restrict
    (
        meshLevel<Type,MeshType>&,
        const meshLevel<Type,MeshType>&,
        const bool = false
    );
};

}

}

}

#ifdef NoRepository
    #include "harmonicVolumeWeightedRestrictionScheme.C"
#endif

#endif
