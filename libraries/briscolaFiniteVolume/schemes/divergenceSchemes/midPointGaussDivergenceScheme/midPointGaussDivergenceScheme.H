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
    public divergenceScheme<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("midPointGauss");

    // Construct from dictionary and mesh
    midPointGaussDivergenceScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    midPointGaussDivergenceScheme(const fvMesh&);

    // Divergence with a flux defined at faces
    virtual tmp<linearSystem<stencil,Type,MeshType>> div
    (
        const meshField<hexScalar,MeshType>&,
        meshField<Type,MeshType>&
    );
};

}

}

}

#ifdef NoRepository
    #include "midPointGaussDivergenceScheme.C"
#endif

#endif
