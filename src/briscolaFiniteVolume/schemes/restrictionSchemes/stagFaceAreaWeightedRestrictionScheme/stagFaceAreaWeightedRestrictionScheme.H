#ifndef stagFaceAreaWeightedRestrictionScheme_H
#define stagFaceAreaWeightedRestrictionScheme_H

#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class stagFaceAreaWeightedRestrictionScheme
:
    public restrictionScheme<FaceSpace<Type>,staggered>
{
public:

    // Runtime type information
    TypeName("faceAreaWeighted");

    // Construct from stream and mesh
    stagFaceAreaWeightedRestrictionScheme
    (
        const fvMesh&,
        Istream&
    );

    // In-place restriction from a fine cell level to an adjacent coarse cell
    // level
    virtual void restrict
    (
        meshDirection<FaceSpace<Type>,staggered>&,
        const meshDirection<FaceSpace<Type>,staggered>&,
        const bool = false
    );
};

}

}

}

#ifdef NoRepository
    #include "stagFaceAreaWeightedRestrictionScheme.C"
#endif

#endif
