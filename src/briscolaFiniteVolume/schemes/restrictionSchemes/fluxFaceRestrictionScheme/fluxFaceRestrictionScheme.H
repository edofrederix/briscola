#ifndef fluxFaceRestrictionScheme_H
#define fluxFaceRestrictionScheme_H

#include "faceRestrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Forward declaration

template<class Type, class MeshType>
class fluxFaceRestrictionScheme;


// Colocated

template<class Type>
class fluxFaceRestrictionScheme<Type,colocated>
:
    public faceRestrictionScheme<Type,colocated>
{
public:

    // Runtime type information
    TypeName("flux");

    // Use constructors from base
    using faceRestrictionScheme<Type,colocated>::faceRestrictionScheme;

    // In-place restriction from a fine face level to an adjacent coarse face
    // level
    virtual void restrict
    (
        faceField<Type,colocated>&,
        const label l
    );
};


// Staggered

template<class Type>
class fluxFaceRestrictionScheme<Type,staggered>
:
    public faceRestrictionScheme<Type,staggered>
{
public:

    // Runtime type information
    TypeName("flux");

    // Use constructors from base
    using faceRestrictionScheme<Type,staggered>::faceRestrictionScheme;

    // In-place restriction from a fine face level to an adjacent coarse face
    // level
    virtual void restrict
    (
        faceField<Type,staggered>&,
        const label l
    );
};

}

}

}

#ifdef NoRepository
    #include "fluxFaceRestrictionScheme.C"
#endif

#endif
