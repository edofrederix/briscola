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

    // Construct from stream and mesh
    linearGaussGradientScheme(const fvMesh&, Istream& = scheme::nullStream);

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
};

}

}

}

#ifdef NoRepository
    #include "linearGaussGradientScheme.C"
#endif

#endif
