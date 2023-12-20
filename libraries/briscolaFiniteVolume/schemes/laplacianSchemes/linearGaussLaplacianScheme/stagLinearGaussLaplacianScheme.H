#ifndef stagLinearGaussLaplacianScheme_H
#define stagLinearGaussLaplacianScheme_H

#include "linearGaussLaplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class stagLinearGaussLaplacianScheme
:
    public linearGaussLaplacianScheme<stencil,Type,staggered>
{
public:

    // Construct from dictionary and mesh
    stagLinearGaussLaplacianScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    stagLinearGaussLaplacianScheme(const fvMesh&);

    // Implicit Laplacian
    virtual tmp<linearSystem<stencil,Type,staggered>> imLaplacian
    (
        const meshField<lowerFaceScalar,staggered>* lambdaPtr,
        meshField<Type,staggered>& field
    );
};

}

}

}

#ifdef NoRepository
    #include "stagLinearGaussLaplacianScheme.C"
#endif

#endif
