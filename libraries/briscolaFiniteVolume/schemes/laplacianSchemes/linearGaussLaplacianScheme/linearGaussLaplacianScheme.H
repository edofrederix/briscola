#ifndef linearGaussLaplacianScheme_H
#define linearGaussLaplacianScheme_H

#include "laplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class linearGaussLaplacianScheme
:
    public laplacianScheme<SType,Type,MeshType>
{
public:

    // Runtime type information
    TypeName("linearGauss");

    // Construct from dictionary and mesh
    linearGaussLaplacianScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    linearGaussLaplacianScheme(const fvMesh&);

    // Explicit laplacian schemes
    virtual tmp<meshField<Type,MeshType>> exLaplacian
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
    #include "linearGaussLaplacianScheme.C"
#endif

#endif
