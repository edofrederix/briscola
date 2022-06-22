#ifndef EulerDdtScheme_H
#define EulerDdtScheme_H

#include "ddtScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class EulerDdtScheme
:
    public ddtScheme<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("Euler");

    // Construct from dictionary and mesh
    EulerDdtScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    EulerDdtScheme(const fvMesh&);

    // Linear system from ddt(u)
    virtual tmp<linearSystem<diagStencil,Type,MeshType>> ddt
    (
        meshField<Type,MeshType>&
    );

    // Linear system from ddt(coeff,u)
    virtual tmp<linearSystem<diagStencil,Type,MeshType>> ddt
    (
        const meshField<scalar,MeshType>&,
        meshField<Type,MeshType>&
    );
};

}

}

}

#ifdef NoRepository
    #include "EulerDdtScheme.C"
#endif

#endif
