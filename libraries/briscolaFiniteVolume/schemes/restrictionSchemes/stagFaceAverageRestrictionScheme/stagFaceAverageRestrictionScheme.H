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
    public restrictionScheme<Type,staggered>
{
public:

    // Runtime type information
    TypeName("faceAverage");

    // Construct from dictionary and mesh
    stagFaceAverageRestrictionScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    stagFaceAverageRestrictionScheme(const fvMesh&);

    // In-place restriction from a fine cell level to an adjacent coarse cell
    // level
    virtual void restrict
    (
        meshDirection<Type,staggered>&,
        const meshDirection<Type,staggered>&,
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
