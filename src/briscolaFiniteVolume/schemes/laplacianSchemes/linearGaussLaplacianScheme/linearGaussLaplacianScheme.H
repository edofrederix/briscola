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
    public laplacianScheme<stencil,Type,MeshType>
{
public:

    // Runtime type information
    TypeName("linearGauss");

    // Construct from stream and mesh
    linearGaussLaplacianScheme(const fvMesh&, Istream& = scheme::nullStream);

    // Implicit laplacian scheme
    virtual tmp<linearSystem<stencil,Type,MeshType>> imLaplacian
    (
        const faceField<scalar,MeshType>* lambdaPtr,
        const meshField<Type,MeshType>& field
    );

    // Explicit laplacian schemes
    virtual tmp<meshField<Type,MeshType>> exLaplacian
    (
        const faceField<scalar,MeshType>* lambdaPtr,
        const meshField<Type,MeshType>& field
    );
};

}

}

}

#ifdef NoRepository
    #include "linearGaussLaplacianScheme.C"
#endif

#endif
