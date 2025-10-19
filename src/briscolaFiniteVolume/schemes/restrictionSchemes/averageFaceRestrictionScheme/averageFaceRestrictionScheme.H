#ifndef averageFaceRestrictionScheme_H
#define averageFaceRestrictionScheme_H

#include "faceRestrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Forward declaration

template<class Type, class MeshType>
class averageFaceRestrictionScheme;


// Colocated

template<class Type>
class averageFaceRestrictionScheme<Type,colocated>
:
    public faceRestrictionScheme<Type,colocated>
{
public:

    // Runtime type information
    TypeName("average");

    // Use constructors from base
    using faceRestrictionScheme<Type,colocated>::faceRestrictionScheme;

    // In-place restriction from a fine face level to an adjacent coarse face
    // level
    virtual void restrict
    (
        faceField<Type,colocated>&,
        const label l,
        const label d
    );
};


// Staggered

template<class Type>
class averageFaceRestrictionScheme<Type,staggered>
:
    public faceRestrictionScheme<Type,staggered>
{
public:

    // Runtime type information
    TypeName("average");

    // Use constructors from base
    using faceRestrictionScheme<Type,staggered>::faceRestrictionScheme;

    // In-place restriction from a fine face level to an adjacent coarse face
    // level
    virtual void restrict
    (
        faceField<Type,staggered>&,
        const label l,
        const label d
    );
};

}

}

}

#ifdef NoRepository
    #include "averageFaceRestrictionScheme.C"
#endif

#endif
