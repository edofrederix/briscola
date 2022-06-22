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
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& fvMsh
        ),
        (dict, fvMsh)
    );

    // Construct from dictionary and mesh
    ddtScheme(const dictionary&, const fvMesh&);

    // Copy constructor
    ddtScheme(const ddtScheme<Type,MeshType>&);

    // Selector
    static autoPtr<ddtScheme<Type,MeshType>> New(const word, const fvMesh&);

    // Destructor
    virtual ~ddtScheme();

    // Time step size
    inline scalar deltaT() const { return this->fvMsh().time().deltaTValue(); }

    // Linear system from ddt(u)
    virtual tmp<linearSystem<diagStencil,Type,MeshType>> ddt
    (
        meshField<Type,MeshType>&
    ) = 0;

    // Linear system from ddt(coeff,u)
    virtual tmp<linearSystem<diagStencil,Type,MeshType>> ddt
    (
        const meshField<scalar,MeshType>&,
        meshField<Type,MeshType>&
    ) = 0;
};

}

}

}

#ifdef NoRepository
    #include "ddtScheme.C"
#endif

#endif
