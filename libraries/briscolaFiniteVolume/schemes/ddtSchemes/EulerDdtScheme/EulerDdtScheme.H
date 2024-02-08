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

    // Pure virtual functions, to be defined by derived schemes

    virtual tmp<linearSystem<diagStencil,Type,MeshType>> imDdt
    (
        const meshField<scalar,MeshType>* lambdaPtr,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    );

    virtual tmp<meshField<Type,MeshType>> exDdt
    (
        const meshField<scalar,MeshType>* lambdaPtr,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    );
};

}

}

}

#ifdef NoRepository
    #include "EulerDdtScheme.C"
#endif

#endif
