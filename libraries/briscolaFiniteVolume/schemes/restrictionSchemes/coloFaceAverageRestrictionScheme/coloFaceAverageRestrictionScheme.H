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

    // Construct from stream and mesh
    coloFaceAverageRestrictionScheme
    (
        const fvMesh&,
        Istream& = scheme::nullStream
    );

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
