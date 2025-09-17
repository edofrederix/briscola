#ifndef averageRestrictionScheme_H
#define averageRestrictionScheme_H

#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class averageRestrictionScheme
:
    public restrictionScheme<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("average");

    // Use base constructors
    using restrictionScheme<Type,MeshType>::restrictionScheme;

    // In-place restriction from a fine cell level to an adjacent coarse cell
    // level, with option to scale by cell volume
    virtual void restrict
    (
        meshDirection<Type,MeshType>&,
        const meshDirection<Type,MeshType>&,
        const bool = false
    );
};

}

}

}

#ifdef NoRepository
    #include "averageRestrictionScheme.C"
#endif

#endif
