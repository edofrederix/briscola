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

template<class Type, class MeshType>
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
    divergenceScheme(const divergenceScheme<Type,MeshType>&);

    // Selector
    static autoPtr<divergenceScheme<Type,MeshType>> New(const word, const fvMesh&);

    // Destructor
    virtual ~divergenceScheme();

    // Divergence with a flux defined at faces
    virtual tmp<linearSystem<stencil,Type,MeshType>> div
    (
        const meshField<faceScalar,MeshType>&,
        meshField<Type,MeshType>&
    ) = 0;
};

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>> explicitDiv
(
    const meshField<FaceSpace<Type>,MeshType>&
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
