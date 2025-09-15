#ifndef coloFaceAreaWeightedRestrictionScheme_H
#define coloFaceAreaWeightedRestrictionScheme_H

#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class coloFaceAreaWeightedRestrictionScheme
:
    public restrictionScheme<FaceSpace<Type>,colocated>
{
public:

    // Runtime type information
    TypeName("faceAreaWeighted");

    // Construct from stream and mesh
    coloFaceAreaWeightedRestrictionScheme
    (
        const fvMesh&,
        Istream& = scheme::nullStream
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
    #include "coloFaceAreaWeightedRestrictionScheme.C"
#endif

#endif
