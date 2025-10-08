#ifndef faceAreaWeightedRestrictionScheme_H
#define faceAreaWeightedRestrictionScheme_H

#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Forward declaration

template<class Type, class MeshType>
class faceAreaWeightedRestrictionScheme;


// Colocated

template<class Type>
class faceAreaWeightedRestrictionScheme<Type,colocated>
:
    public restrictionScheme<FaceSpace<Type>,colocated>
{
public:

    // Runtime type information
    TypeName("faceAreaWeighted");

    // Use constructors from base
    using restrictionScheme<FaceSpace<Type>,colocated>::restrictionScheme;

    // In-place restriction from a fine cell level to an adjacent coarse cell
    // level
    virtual void restrict
    (
        meshDirection<FaceSpace<Type>,colocated>&,
        const meshDirection<FaceSpace<Type>,colocated>&,
        const bool = false
    );
};


// Staggered

template<class Type>
class faceAreaWeightedRestrictionScheme<Type,staggered>
:
    public restrictionScheme<FaceSpace<Type>,staggered>
{
public:

    // Runtime type information
    TypeName("faceAreaWeighted");

    // Use constructors from base
    using restrictionScheme<FaceSpace<Type>,staggered>::restrictionScheme;

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
    #include "faceAreaWeightedRestrictionScheme.C"
#endif

#endif
