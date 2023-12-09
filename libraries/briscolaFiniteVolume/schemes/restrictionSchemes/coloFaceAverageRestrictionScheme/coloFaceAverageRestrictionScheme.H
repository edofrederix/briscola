#ifndef coloFaceAverageRestrictionScheme_H
#define coloFaceAverageRestrictionScheme_H

#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class coloFaceAverageRestrictionScheme
:
    public restrictionScheme<Type,colocated>
{
public:

    // Runtime type information
    TypeName("faceAverage");

    // Construct from dictionary and mesh
    coloFaceAverageRestrictionScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    coloFaceAverageRestrictionScheme(const fvMesh&);

    // In-place restriction from a fine cell level to an adjacent coarse cell
    // level
    virtual void restrict
    (
        meshDirection<Type,colocated>&,
        const meshDirection<Type,colocated>&,
        const bool = false
    );
};

}

}

}

#ifdef NoRepository
    #include "coloFaceAverageRestrictionScheme.C"
#endif

#endif
