#ifndef harmonicFaceAreaWeightedRestrictionScheme_H
#define harmonicFaceAreaWeightedRestrictionScheme_H

#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Forward declaration

template<class Type, class MeshType>
class harmonicFaceAreaWeightedRestrictionScheme;


// Colocated

template<class Type>
class harmonicFaceAreaWeightedRestrictionScheme<Type,colocated>
:
    public restrictionScheme<FaceSpace<Type>,colocated>
{
public:

    // Runtime type information
    TypeName("harmonicFaceAreaWeighted");

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
class harmonicFaceAreaWeightedRestrictionScheme<Type,staggered>
:
    public restrictionScheme<FaceSpace<Type>,staggered>
{
public:

    // Runtime type information
    TypeName("harmonicFaceAreaWeighted");

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
    #include "harmonicFaceAreaWeightedRestrictionScheme.C"
#endif

#endif
