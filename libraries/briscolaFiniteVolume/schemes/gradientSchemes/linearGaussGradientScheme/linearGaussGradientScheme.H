#ifndef linearGaussGradientScheme_H
#define linearGaussGradientScheme_H

#include "gradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class linearGaussGradientScheme
:
    public gradientScheme<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("linearGauss");

    // Construct from dictionary and mesh
    linearGaussGradientScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    linearGaussGradientScheme(const fvMesh&);

    // Gradient of u from same mesh type
    virtual tmp
    <
        meshField
        <
            typename gradientScheme<Type,MeshType>::GradType,
            MeshType
        >
    >
    grad
    (
        const meshField<Type,MeshType>&
    );

    // Staggered gradient of a colocated field, where each mesh direction holds
    // a component of the gradient
    virtual tmp<meshField<Type,staggered>> stagGrad
    (
        const meshField<Type,colocated>&
    );
};

}

}

}

#ifdef NoRepository
    #include "linearGaussGradientScheme.C"
#endif

#endif
