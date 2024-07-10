#ifndef gradientScheme_H
#define gradientScheme_H

#include "scheme.H"
#include "meshFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class gradientScheme
:
    public scheme
{
public:

    // Type of the gradient
    typedef typename outerProduct<vector, Type>::type GradType;

    // Runtime type information
    TypeName("gradientScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        gradientScheme,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& fvMsh
        ),
        (dict, fvMsh)
    );

    // Construct from dictionary and mesh
    gradientScheme(const dictionary&, const fvMesh&);

    // Copy constructor
    gradientScheme(const gradientScheme<Type,MeshType>&);

    // Selector
    static autoPtr<gradientScheme<Type,MeshType>> New(const word, const fvMesh&);

    // Destructor
    virtual ~gradientScheme();

    // Gradient of u from same mesh type
    virtual tmp<meshField<GradType,MeshType>> grad
    (
        const meshField<Type,MeshType>&
    ) = 0;

    // Staggered gradient of a colocated field, where each mesh direction in the
    // staggered field holds a component of the gradient
    virtual tmp<meshField<Type,staggered>> stagGrad
    (
        const meshField<Type,colocated>&
    ) = 0;
};

// Gradient of a face field

template<class Type, class MeshType>
tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
grad
(
    const meshField<FaceSpace<Type>,MeshType>&
);

template<class Type, class MeshType>
tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
grad
(
    const tmp<meshField<FaceSpace<Type>,MeshType>>&
);

template<class Type, class MeshType>
tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
grad
(
    const meshField<LowerFaceSpace<Type>,MeshType>&
);

template<class Type, class MeshType>
tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
grad
(
    const tmp<meshField<LowerFaceSpace<Type>,MeshType>>&
);

}

}

}

#ifdef NoRepository
    #include "gradientScheme.C"
#endif

#endif
