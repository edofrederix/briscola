#ifndef midPointGaussGradientScheme_H
#define midPointGaussGradientScheme_H

#include "gradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class midPointGaussGradientScheme
:
    public gradientScheme<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("midPointGauss");

    // Construct from stream and mesh
    midPointGaussGradientScheme(const fvMesh&, Istream& = scheme::nullStream);

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
    #include "midPointGaussGradientScheme.C"
#endif

#endif
