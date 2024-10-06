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

    // Construct from stream and mesh
    harmonicVolumeWeightedRestrictionScheme
    (
        const fvMesh&,
        Istream& = scheme::nullStream
    );

    // In-place restriction from a fine cell level to an adjacent coarse cell
    // level, with option to scale by cell volume
    virtual void restrict
    (
        meshDirection<Type,MeshType>&,
        const meshDirection<Type,MeshType>&,
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
