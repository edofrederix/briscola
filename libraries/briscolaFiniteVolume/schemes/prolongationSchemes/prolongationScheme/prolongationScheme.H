#ifndef prolongationScheme_H
#define prolongationScheme_H

#include "scheme.H"
#include "meshFields.H"
#include "ops.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class prolongationScheme
:
    public scheme
{
public:

    // Runtime type information
    TypeName("prolongationScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        prolongationScheme,
        Istream,
        (
            const fvMesh& fvMsh,
            Istream& is
        ),
        (fvMsh, is)
    );

    // Construct from stream and mesh
    prolongationScheme(const fvMesh&, Istream& = scheme::nullStream);

    // Copy constructor
    prolongationScheme(const prolongationScheme<Type,MeshType>&);

    // Selector
    static autoPtr<prolongationScheme<Type,MeshType>> New
    (
        const fvMesh&,
        const word
    );

    // Destructor
    virtual ~prolongationScheme();

    // In-place prolongation from a coarse cell level to the adjacent fine cell
    // level, using the = operation
    virtual void prolong
    (
        meshDirection<Type,MeshType>&,
        const meshDirection<Type,MeshType>&,
        const eqOp<Type>&
    ) = 0;

    // In-place prolongation from a coarse cell level to the adjacent fine cell
    // level, using the += operation
    virtual void prolong
    (
        meshDirection<Type,MeshType>&,
        const meshDirection<Type,MeshType>&,
        const plusEqOp<Type>&
    ) = 0;

    // In-place prolongation from a coarse cell level to the adjacent fine cell
    // level, using the -= operation
    virtual void prolong
    (
        meshDirection<Type,MeshType>&,
        const meshDirection<Type,MeshType>&,
        const minusEqOp<Type>&
    ) = 0;
};

}

}

}

#ifdef NoRepository
    #include "prolongationScheme.C"
#endif

#endif
