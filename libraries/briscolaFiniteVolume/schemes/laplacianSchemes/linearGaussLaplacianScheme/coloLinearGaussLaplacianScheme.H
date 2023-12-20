#ifndef coloLinearGaussLaplacianScheme_H
#define coloLinearGaussLaplacianScheme_H

#include "linearGaussLaplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class coloLinearGaussLaplacianScheme
:
    public linearGaussLaplacianScheme<symmStencil,Type,colocated>
{
public:

    // Construct from dictionary and mesh
    coloLinearGaussLaplacianScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    coloLinearGaussLaplacianScheme(const fvMesh&);

    // Implicit Laplacian
    virtual tmp<linearSystem<symmStencil,Type,colocated>> imLaplacian
    (
        const meshField<lowerFaceScalar,colocated>* lambdaPtr,
        meshField<Type,colocated>& field
    );
};

}

}

}

#ifdef NoRepository
    #include "coloLinearGaussLaplacianScheme.C"
#endif

#endif
