#ifndef fluxRestrictionScheme_H
#define fluxRestrictionScheme_H

#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Forward declaration

template<class Type, class MeshType>
class fluxRestrictionScheme;


// Colocated

template<class Type>
class fluxRestrictionScheme<Type,colocated>
:
    public restrictionScheme<FaceSpace<Type>,colocated>
{
public:

    // Runtime type information
    TypeName("flux");

    // Use constructors from base
    using restrictionScheme<FaceSpace<Type>,colocated>:: restrictionScheme;

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
    #include "fluxRestrictionScheme.C"
#endif

#endif
