#ifndef linearGaussDivergenceScheme_H
#define linearGaussDivergenceScheme_H

#include "divergenceScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class linearGaussDivergenceScheme
:
    public divergenceScheme<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("linearGauss");

    // Construct from dictionary and mesh
    linearGaussDivergenceScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    linearGaussDivergenceScheme(const fvMesh&);

    // Implicit divergence operator
    virtual tmp<linearSystem<stencil,Type,MeshType>> imDiv
    (
        const meshField<lowerFaceScalar,MeshType>&,
        meshField<Type,MeshType>&
    );

    // Explicit divergence operator
    virtual tmp<meshField<Type,MeshType>> exDiv
    (
        const meshField<lowerFaceScalar,MeshType>&,
        meshField<Type,MeshType>&
    );
};

}

}

}

#ifdef NoRepository
    #include "linearGaussDivergenceScheme.C"
#endif

#endif
