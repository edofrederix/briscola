#ifndef fluxRestrictionScheme_H
#define fluxRestrictionScheme_H

#include "restrictionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
class fluxRestrictionScheme
:
    public restrictionScheme<faceScalar,MeshType>
{
public:

    // Runtime type information
    TypeName("flux");

    // Construct from dictionary and mesh
    fluxRestrictionScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    fluxRestrictionScheme(const fvMesh&);

    // In-place restriction from a fine cell level to an adjacent coarse cell
    // level
    virtual void restrict
    (
        meshDirection<faceScalar,MeshType>&,
        const meshDirection<faceScalar,MeshType>&,
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
