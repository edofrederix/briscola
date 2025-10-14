#ifndef midPointGaussDivergenceScheme_H
#define midPointGaussDivergenceScheme_H

#include "divergenceScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class midPointGaussDivergenceScheme
:
    public divergenceScheme<stencil,Type,MeshType>
{
public:

    // Runtime type information
    TypeName("midPointGauss");

    // Construct from stream and mesh
    midPointGaussDivergenceScheme(const fvMesh&, Istream& = scheme::nullStream);

    // Implicit divergence scheme
    virtual tmp<linearSystem<stencil,Type,MeshType>> imDiv
    (
        const meshField<faceScalar,MeshType>& phi,
        const meshField<Type,MeshType>& field
    );

    // Explicit divergence scheme
    virtual tmp<meshField<Type,MeshType>> exDiv
    (
        const meshField<faceScalar,MeshType>& phi,
        const meshField<Type,MeshType>& field
    );
};

}

}

}

#ifdef NoRepository
    #include "midPointGaussDivergenceScheme.C"
#endif

#endif
