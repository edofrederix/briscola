#ifndef linearGaussLaplacianScheme_H
#define linearGaussLaplacianScheme_H

#include "laplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class linearGaussLaplacianScheme
:
    public laplacianScheme<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("linearGauss");

    // Construct from dictionary and mesh
    linearGaussLaplacianScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    linearGaussLaplacianScheme(const fvMesh&);

    // Implicit Laplacian
    virtual tmp<linearSystem<stencil,Type,MeshType>> imLaplacian
    (
        const meshField<lowerFaceScalar,MeshType>*,
        meshField<Type,MeshType>&
    );

    // Explicit Laplacian
    virtual tmp<meshField<Type,MeshType>> exLaplacian
    (
        const meshField<lowerFaceScalar,MeshType>*,
        meshField<Type,MeshType>&
    );
};

}

}

}

#ifdef NoRepository
    #include "linearGaussLaplacianScheme.C"
#endif

#endif
