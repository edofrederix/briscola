#ifndef laplacianScheme_H
#define laplacianScheme_H

#include "scheme.H"
#include "linearSystem.H"
#include "stencil.H"
#include "ops.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
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
        Istream,
        (
            const fvMesh& fvMsh,
            Istream& is
        ),
        (fvMsh, is)
    );

    // Construct from stream and mesh
    laplacianScheme(const fvMesh&, Istream& = scheme::nullStream);

    // Copy constructor
    laplacianScheme(const laplacianScheme<SType,Type,MeshType>&);

    // Selector
    static autoPtr<laplacianScheme<SType,Type,MeshType>> New
    (
        const fvMesh&,
        const word
    );

    // Destructor
    virtual ~laplacianScheme();

    // Pure virtual functions, to be defined by derived schemes

        virtual tmp<linearSystem<SType,Type,MeshType>> imLaplacian
        (
            const meshField<faceScalar,MeshType>* lambdaPtr,
            const meshField<Type,MeshType>& field
        ) = 0;

        virtual tmp<meshField<Type,MeshType>> exLaplacian
        (
            const meshField<faceScalar,MeshType>* lambdaPtr,
            const meshField<Type,MeshType>& field
        ) = 0;

    // Support functions relying on the pure virtual functions above

        // Implicit

            // Return

                // Coefficient as reference
                inline tmp<linearSystem<SType,Type,MeshType>> imLaplacian
                (
                    const meshField<faceScalar,MeshType>& lambda,
                    const meshField<Type,MeshType>& field
                );

                // Coefficient as tmp
                inline tmp<linearSystem<SType,Type,MeshType>> imLaplacian
                (
                    const tmp<meshField<faceScalar,MeshType>>& tLambda,
                    const meshField<Type,MeshType>& field
                );

                // Constant coefficient
                inline tmp<linearSystem<SType,Type,MeshType>> imLaplacian
                (
                    const scalar lambda,
                    const meshField<Type,MeshType>& field
                );

                // No coefficient
                inline tmp<linearSystem<SType,Type,MeshType>> imLaplacian
                (
                    const meshField<Type,MeshType>& field
                );

        // Explicit

            // Return

                // Coefficient as reference
                inline tmp<meshField<Type,MeshType>> exLaplacian
                (
                    const meshField<faceScalar,MeshType>& lambda,
                    const meshField<Type,MeshType>& field
                );

                // Coefficient as tmp
                inline tmp<meshField<Type,MeshType>> exLaplacian
                (
                    const tmp<meshField<faceScalar,MeshType>>& tLambda,
                    const meshField<Type,MeshType>& field
                );

                // Constant coefficient
                inline tmp<meshField<Type,MeshType>> exLaplacian
                (
                    const scalar lambda,
                    const meshField<Type,MeshType>& field
                );

                // No coefficient
                inline tmp<meshField<Type,MeshType>> exLaplacian
                (
                    const meshField<Type,MeshType>& field
                );
};

}

}

}

#include "laplacianSchemeI.H"

#ifdef NoRepository
    #include "laplacianScheme.C"
#endif

#endif
