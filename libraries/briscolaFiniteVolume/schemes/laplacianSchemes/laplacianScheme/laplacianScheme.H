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

    // Laplacian with a scalar diffusivity defined at faces
    virtual tmp<linearSystem<stencil,Type,MeshType>> laplacian
    (
        const meshField<faceScalar,MeshType>&,
        meshField<Type,MeshType>&
    ) = 0;

    // Laplacian
    virtual tmp<linearSystem<stencil,Type,MeshType>> laplacian
    (
        meshField<Type,MeshType>& field
    )
    {
        meshField<faceScalar,MeshType> lambda("lambda", field.fvMsh());
        lambda = faceScalar::one;
        return this->laplacian(lambda, field);
    }
};

}

}

}

#ifdef NoRepository
    #include "laplacianScheme.C"
#endif

#endif
