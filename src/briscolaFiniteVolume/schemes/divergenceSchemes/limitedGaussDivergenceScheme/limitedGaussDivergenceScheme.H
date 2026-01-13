#ifndef limitedGaussDivergenceScheme_H
#define limitedGaussDivergenceScheme_H

#include "divergenceScheme.H"
#include "limiterScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class limitedGaussDivergenceScheme
:
    public divergenceScheme<stencil,Type,MeshType>
{
    // Limiter scheme
    autoPtr<limiterScheme<Type,MeshType>> limiter_;

public:

    // Runtime type information
    TypeName("limitedGauss");

    // Construct from stream and mesh
    limitedGaussDivergenceScheme(const fvMesh&, Istream& = scheme::nullStream);

    // Implicit divergence scheme
    virtual tmp<linearSystem<stencil,Type,MeshType>> imDiv
    (
        const faceField<scalar,MeshType>& phi,
        const meshField<Type,MeshType>& field
    );

    // Explicit divergence scheme
    virtual tmp<meshField<Type,MeshType>> exDiv
    (
        const faceField<scalar,MeshType>& phi,
        const meshField<Type,MeshType>& field
    );
};

}

}

}

#ifdef NoRepository
    #include "limitedGaussDivergenceScheme.C"
#endif

#endif
