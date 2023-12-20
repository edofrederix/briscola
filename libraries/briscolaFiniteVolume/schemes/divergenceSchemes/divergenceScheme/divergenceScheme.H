#ifndef divergenceScheme_H
#define divergenceScheme_H

#include "scheme.H"
#include "linearSystem.H"
#include "stencil.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class divergenceScheme
:
    public scheme
{
public:

    // Runtime type information
    TypeName("divergenceScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        divergenceScheme,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& fvMsh
        ),
        (dict, fvMsh)
    );

    // Construct from dictionary and mesh
    divergenceScheme(const dictionary&, const fvMesh&);

    // Copy constructor
    divergenceScheme(const divergenceScheme<SType,Type,MeshType>&);

    // Selector
    static autoPtr<divergenceScheme<SType,Type,MeshType>> New
    (
        const word,
        const fvMesh&
    );

    // Destructor
    virtual ~divergenceScheme();

    // Implicit divergence operator
    virtual tmp<linearSystem<SType,Type,MeshType>> imDiv
    (
        const meshField<lowerFaceScalar,MeshType>&,
        meshField<Type,MeshType>&
    ) = 0;

    // Explicit divergence operator
    virtual tmp<meshField<Type,MeshType>> exDiv
    (
        const meshField<lowerFaceScalar,MeshType>&,
        meshField<Type,MeshType>&
    ) = 0;
};

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>> explicitDiv
(
    const meshField<LowerFaceSpace<Type>,MeshType>&
);

template<class Type>
tmp<meshField<Type,colocated>> explicitColoDiv
(
    const meshField<Type,staggered>&
);

}

}

}

#ifdef NoRepository
    #include "divergenceScheme.C"
#endif

#endif
