#ifndef stencilLinearGaussLaplacianScheme_H
#define stencilLinearGaussLaplacianScheme_H

#include "linearGaussLaplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class stencilLinearGaussLaplacianScheme
:
    public linearGaussLaplacianScheme<stencil,Type,MeshType>
{
public:

    // Construct from dictionary and mesh
    stencilLinearGaussLaplacianScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    stencilLinearGaussLaplacianScheme(const fvMesh&);

    // Implicit laplacian scheme
    virtual tmp<linearSystem<stencil,Type,MeshType>> imLaplacian
    (
        const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    );
};

}

}

}

#ifdef NoRepository
    #include "stencilLinearGaussLaplacianScheme.C"
#endif

#endif
