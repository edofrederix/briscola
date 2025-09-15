#ifndef stagHarmonicFaceAreaWeightedRestrictionScheme_H
#define stagHarmonicFaceAreaWeightedRestrictionScheme_H

#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class stagHarmonicFaceAreaWeightedRestrictionScheme
:
    public restrictionScheme<FaceSpace<Type>,staggered>
{
public:

    // Runtime type information
    TypeName("harmonicFaceAreaWeighted");

    // Construct from stream and mesh
    stagHarmonicFaceAreaWeightedRestrictionScheme
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
    #include "stagHarmonicFaceAreaWeightedRestrictionScheme.C"
#endif

#endif
