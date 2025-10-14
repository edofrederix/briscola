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

    // Construct from stream and mesh
    linearGaussLaplacianScheme(const fvMesh&, Istream& = scheme::nullStream);

    // Explicit laplacian schemes
    virtual tmp<meshField<Type,MeshType>> exLaplacian
    (
        const meshField<faceScalar,MeshType>* lambdaPtr,
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
