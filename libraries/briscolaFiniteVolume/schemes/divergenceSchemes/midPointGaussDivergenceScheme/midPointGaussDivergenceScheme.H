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

    // Construct from dictionary and mesh
    midPointGaussDivergenceScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    midPointGaussDivergenceScheme(const fvMesh&);

    // Implicit divergence scheme
    virtual tmp<linearSystem<stencil,Type,MeshType>> imDiv
    (
        const meshField<lowerFaceScalar,MeshType>& phi,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    );

    // Explicit divergence scheme
    virtual tmp<meshField<Type,MeshType>> exDiv
    (
        const meshField<lowerFaceScalar,MeshType>& phi,
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
