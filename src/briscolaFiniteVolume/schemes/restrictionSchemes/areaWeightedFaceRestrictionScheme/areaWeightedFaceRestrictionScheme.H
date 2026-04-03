#ifndef areaWeightedFaceRestrictionScheme_H
#define areaWeightedFaceRestrictionScheme_H

#include "faceRestrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Forward declaration

template<class Type, class MeshType>
class areaWeightedFaceRestrictionScheme;


// Colocated

template<class Type>
class areaWeightedFaceRestrictionScheme<Type,colocated>
:
    public faceRestrictionScheme<Type,colocated>
{
public:

    // Runtime type information
    TypeName("areaWeighted");

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
class areaWeightedFaceRestrictionScheme<Type,staggered>
:
    public faceRestrictionScheme<Type,staggered>
{
public:

    // Runtime type information
    TypeName("areaWeighted");

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
    #include "areaWeightedFaceRestrictionScheme.C"
#endif

#endif
