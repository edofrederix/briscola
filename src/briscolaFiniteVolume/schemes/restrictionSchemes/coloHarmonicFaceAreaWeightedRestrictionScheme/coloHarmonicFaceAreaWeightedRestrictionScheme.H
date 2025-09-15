#ifndef coloHarmonicFaceAreaWeightedRestrictionScheme_H
#define coloHarmonicFaceAreaWeightedRestrictionScheme_H

#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class coloHarmonicFaceAreaWeightedRestrictionScheme
:
    public restrictionScheme<FaceSpace<Type>,colocated>
{
public:

    // Runtime type information
    TypeName("harmonicFaceAreaWeighted");

    // Construct from stream and mesh
    coloHarmonicFaceAreaWeightedRestrictionScheme
    (
        const fvMesh&,
        Istream&
    );

    // In-place restriction from a fine cell level to an adjacent coarse cell
    // level
    virtual void restrict
    (
        meshDirection<FaceSpace<Type>,colocated>&,
        const meshDirection<FaceSpace<Type>,colocated>&,
        const bool = false
    );
};

}

}

}

#ifdef NoRepository
    #include "coloHarmonicFaceAreaWeightedRestrictionScheme.C"
#endif

#endif
