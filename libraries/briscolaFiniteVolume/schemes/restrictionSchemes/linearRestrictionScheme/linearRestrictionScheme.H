#ifndef linearRestrictionScheme_H
#define linearRestrictionScheme_H

#include "restrictionScheme.H"
#include "PtrList.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class linearRestrictionScheme
:
    public restrictionScheme<Type,MeshType>
{
    // Cached weights
    meshField<vertexScalar,MeshType> weights_;

    // Set weights
    void setWeights();

public:

    // Runtime type information
    TypeName("linear");

    // Construct from dictionary and mesh
    linearRestrictionScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    linearRestrictionScheme(const fvMesh&);

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
    #include "linearRestrictionScheme.C"
#endif

#endif
