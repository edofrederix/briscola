#ifndef upwindDivergenceScheme_H
#define upwindDivergenceScheme_H

#include "divergenceScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class upwindDivergenceScheme
:
    public divergenceScheme<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("upwind");

    // Construct from dictionary and mesh
    upwindDivergenceScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    upwindDivergenceScheme(const fvMesh&);

    // Divergence with a flux defined at faces
    virtual tmp<linearSystem<stencil,Type,MeshType>> div
    (
        const meshField<faceScalar,MeshType>&,
        meshField<Type,MeshType>&
    );
};

}

}

}

#ifdef NoRepository
    #include "upwindDivergenceScheme.C"
#endif

#endif
