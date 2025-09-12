#ifndef stagFaceAverageRestrictionScheme_H
#define stagFaceAverageRestrictionScheme_H

#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class stagFaceAverageRestrictionScheme
:
    public restrictionScheme<FaceSpace<Type>,staggered>
{
public:

    // Runtime type information
    TypeName("faceAverage");

    // Construct from stream and mesh
    stagFaceAverageRestrictionScheme
    (
        const fvMesh&,
        Istream& = scheme::nullStream
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
    #include "stagFaceAverageRestrictionScheme.C"
#endif

#endif
