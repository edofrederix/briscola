#ifndef ddtScheme_H
#define ddtScheme_H

#include "scheme.H"
#include "linearSystem.H"
#include "diagStencil.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class ddtScheme
:
    public scheme
{
public:

    // Runtime type information
    TypeName("ddtScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        ddtScheme,
        Istream,
        (
            const fvMesh& fvMsh,
            Istream& is
        ),
        (fvMsh, is)
    );

    // Construct from stream and mesh
    ddtScheme(const fvMesh&, Istream& = scheme::nullStream);

    // Copy constructor
    ddtScheme(const ddtScheme<Type,MeshType>&);

    // Selector
    static autoPtr<ddtScheme<Type,MeshType>> New(const fvMesh&, const word);

    // Destructor
    virtual ~ddtScheme();

    // Time step size
    inline scalar deltaT() const { return this->fvMsh().time().deltaTValue(); }

    // Pure virtual functions, to be defined by derived schemes

    virtual tmp<linearSystem<diagStencil,Type,MeshType>> imDdt
    (
        const meshField<scalar,MeshType>* lambdaPtr,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    ) = 0;

    virtual tmp<meshField<Type,MeshType>> exDdt
    (
        const meshField<scalar,MeshType>* lambdaPtr,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    ) = 0;

    // Support functions relying on the pure virtual functions above

        // Implicit

            // Coefficient as reference
            inline tmp<linearSystem<diagStencil,Type,MeshType>> imDdt
            (
                const meshField<scalar,MeshType>& lambda,
                const meshField<Type,MeshType>& field,
                const scalar factor = 1.0
            );

            // Coefficient as tmp
            inline tmp<linearSystem<diagStencil,Type,MeshType>> imDdt
            (
                const tmp<meshField<scalar,MeshType>>& tLambda,
                const meshField<Type,MeshType>& field,
                const scalar factor = 1.0

            );

            // Constant coefficient
            inline tmp<linearSystem<diagStencil,Type,MeshType>> imDdt
            (
                const scalar lambda,
                const meshField<Type,MeshType>& field,
                const scalar factor = 1.0
            );

            // No coefficient
            inline tmp<linearSystem<diagStencil,Type,MeshType>> imDdt
            (
                const meshField<Type,MeshType>& field,
                const scalar factor = 1.0
            );

        // Explicit

            // Coefficient as reference
            inline tmp<meshField<Type,MeshType>> exDdt
            (
                const meshField<scalar,MeshType>& lambda,
                const meshField<Type,MeshType>& field,
                const scalar factor = 1.0
            );

            // Coefficient as tmp
            inline tmp<meshField<Type,MeshType>> exDdt
            (
                const tmp<meshField<scalar,MeshType>>& tLambda,
                const meshField<Type,MeshType>& field,
                const scalar factor = 1.0
            );

            // Constant coefficient
            inline tmp<meshField<Type,MeshType>> exDdt
            (
                const scalar lambda,
                const meshField<Type,MeshType>& field,
                const scalar factor = 1.0
            );

            // No coefficient
            inline tmp<meshField<Type,MeshType>> exDdt
            (
                const meshField<Type,MeshType>& field,
                const scalar factor = 1.0
            );
};

}

}

}

#include "ddtSchemeI.H"

#ifdef NoRepository
    #include "ddtScheme.C"
#endif

#endif
