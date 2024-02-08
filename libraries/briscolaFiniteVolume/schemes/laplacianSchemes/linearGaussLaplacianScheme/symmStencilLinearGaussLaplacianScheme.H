#ifndef symmStencilLinearGaussLaplacianScheme_H
#define symmStencilLinearGaussLaplacianScheme_H

#include "linearGaussLaplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class symmStencilLinearGaussLaplacianScheme
:
    public linearGaussLaplacianScheme<symmStencil,Type,MeshType>
{
public:

    // Construct from dictionary and mesh
    symmStencilLinearGaussLaplacianScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    symmStencilLinearGaussLaplacianScheme(const fvMesh&);

    // Implicit laplacian scheme
    virtual tmp<linearSystem<symmStencil,Type,MeshType>> imLaplacian
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
    #include "symmStencilLinearGaussLaplacianScheme.C"
#endif

#endif
