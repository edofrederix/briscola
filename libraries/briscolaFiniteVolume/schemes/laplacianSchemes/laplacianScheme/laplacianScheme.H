#ifndef laplacianScheme_H
#define laplacianScheme_H

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
class laplacianScheme
:
    public scheme
{
public:

    // Runtime type information
    TypeName("laplacianScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        laplacianScheme,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& fvMsh
        ),
        (dict, fvMsh)
    );

    // Construct from dictionary and mesh
    laplacianScheme(const dictionary&, const fvMesh&);

    // Copy constructor
    laplacianScheme(const laplacianScheme<Type,MeshType>&);

    // Selector
    static autoPtr<laplacianScheme<Type,MeshType>> New
    (
        const word,
        const fvMesh&
    );

    // Destructor
    virtual ~laplacianScheme();

    // Implicit Laplacian with pointer to coefficient
    virtual tmp<linearSystem<stencil,Type,MeshType>> imLaplacian
    (
        const meshField<lowerFaceScalar,MeshType>*,
        meshField<Type,MeshType>&
    ) = 0;

    // Explicit Laplacian with pointer to coefficient
    virtual tmp<meshField<Type,MeshType>> exLaplacian
    (
        const meshField<lowerFaceScalar,MeshType>*,
        meshField<Type,MeshType>&
    ) = 0;

    // Implicit Laplacian
    inline tmp<linearSystem<stencil,Type,MeshType>> imLaplacian
    (
        const meshField<lowerFaceScalar,MeshType>& lambda,
        meshField<Type,MeshType>& field
    )
    {
        return this->imLaplacian(&lambda, field);
    }

    // Explicit Laplacian
    inline tmp<meshField<Type,MeshType>> exLaplacian
    (
        const meshField<lowerFaceScalar,MeshType>& lambda,
        meshField<Type,MeshType>& field
    )
    {
        return this->exLaplacian(&lambda, field);
    }

    // Implicit Laplacian without coefficient
    inline tmp<linearSystem<stencil,Type,MeshType>> imLaplacian
    (
        meshField<Type,MeshType>& field
    )
    {
        return this->imLaplacian(nullptr, field);
    }

    // Explicit Laplacian without coefficient
    inline tmp<meshField<Type,MeshType>> exLaplacian
    (
        meshField<Type,MeshType>& field
    )
    {
        return this->exLaplacian(nullptr, field);
    }
};

}

}

}

#ifdef NoRepository
    #include "laplacianScheme.C"
#endif

#endif
