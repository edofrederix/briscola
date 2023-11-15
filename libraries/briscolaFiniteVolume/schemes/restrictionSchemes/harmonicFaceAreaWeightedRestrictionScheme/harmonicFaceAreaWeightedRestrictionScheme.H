#ifndef harmonicFaceAreaWeightedRestrictionScheme_H
#define harmonicFaceAreaWeightedRestrictionScheme_H

#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class harmonicFaceAreaWeightedRestrictionScheme
:
    public restrictionScheme<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("harmonicFaceAreaWeighted");

    // Construct from dictionary and mesh
    harmonicFaceAreaWeightedRestrictionScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    harmonicFaceAreaWeightedRestrictionScheme(const fvMesh&);

    // In-place restriction from a fine cell level to an adjacent coarse cell
    // level
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
    #include "harmonicFaceAreaWeightedRestrictionScheme.C"
#endif

#endif
