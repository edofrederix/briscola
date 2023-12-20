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

    // Implicit Laplacian
    virtual tmp<linearSystem<SType,Type,MeshType>> imLaplacian
    (
        const meshField<lowerFaceScalar,MeshType>*,
        meshField<Type,MeshType>&
    ) = 0;

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
