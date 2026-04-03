#ifndef prolongationScheme_H
#define prolongationScheme_H

#include "scheme.H"

#include "ops.H"

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

    // Selector from type name
    static autoPtr<prolongationScheme<Type,MeshType>> New
    (
        const fvMesh&,
        const word
    );

    // Destructor
    virtual ~prolongationScheme();

    // Default prolongation scheme for this mesh and type
    static const word defaultScheme;

    // In-place prolongation from a coarse cell level to the adjacent fine cell
    // level, using the = operation
    virtual void prolong
    (
        meshLevel<Type,MeshType>&,
        const meshLevel<Type,MeshType>&,
        const eqOp<Type>&
    ) = 0;

    // In-place prolongation from a coarse cell level to the adjacent fine cell
    // level, using the += operation
    virtual void prolong
    (
        meshLevel<Type,MeshType>&,
        const meshLevel<Type,MeshType>&,
        const plusEqOp<Type>&
    ) = 0;

    // In-place prolongation from a coarse cell level to the adjacent fine cell
    // level, using the -= operation
    virtual void prolong
    (
        meshLevel<Type,MeshType>&,
        const meshLevel<Type,MeshType>&,
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
