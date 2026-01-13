#ifndef gradientScheme_H
#define gradientScheme_H

#include "scheme.H"
#include "faceFields.H"
#include "meshFields.H"
#include "meshLevels.H"
#include "meshDirections.H"

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
        Istream,
        (
            const fvMesh& fvMsh,
            Istream& is
        ),
        (fvMsh, is)
    );

    // Construct from stream and mesh
    gradientScheme(const fvMesh&, Istream& = scheme::nullStream);

    // Copy constructor
    gradientScheme(const gradientScheme<Type,MeshType>&);

    // Selector
    static autoPtr<gradientScheme<Type,MeshType>> New
    (
        const fvMesh&,
        const word
    );

    // Destructor
    virtual ~gradientScheme();

    // Gradient of u from same mesh type
    virtual tmp<meshField<GradType,MeshType>> grad
    (
        const meshField<Type,MeshType>&
    ) = 0;
};

// Gradient of a face field

template<class Type, class MeshType>
tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
grad
(
    const faceField<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
grad
(
    const tmp<faceField<Type,MeshType>>&
);

// Staggered gradient of a colocated field, where each mesh direction in the
// staggered field holds a component of the gradient

template<class Type, class MeshType>
tmp<meshField<Type,staggered>> stagGrad
(
    const meshField<Type,colocated>&
);

template<class Type, class MeshType>
tmp<meshField<Type,staggered>> stagGrad
(
    const tmp<meshField<Type,colocated>>&
);

}

}

}

#ifdef NoRepository
    #include "gradientScheme.C"
#endif

#endif
