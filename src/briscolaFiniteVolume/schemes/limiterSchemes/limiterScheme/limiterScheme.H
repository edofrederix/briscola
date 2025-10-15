#ifndef limiterScheme_H
#define limiterScheme_H

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
class limiterScheme
:
    public scheme
{
private:

    // Scalar helper function for the ratio of successive gradients
    tmp<meshField<faceScalar,MeshType>> rType
    (
        const meshField<faceScalar,MeshType>& phi,
        const meshField<scalar,MeshType>& field
    );

    // Vector helper function for the ratio of successive gradients
    tmp<meshField<faceScalar,MeshType>> rType
    (
        const meshField<faceScalar,MeshType>& phi,
        const meshField<vector,MeshType>& field
    );

public:

    // Runtime type information
    TypeName("limiterScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        limiterScheme,
        Istream,
        (
            const fvMesh& fvMsh,
            Istream& is
        ),
        (fvMsh, is)
    );

    // Construct from stream and mesh
    limiterScheme(const fvMesh&, Istream&);

    // Copy constructor
    limiterScheme(const limiterScheme<Type,MeshType>&);

    // Selector from scheme name
    static autoPtr<limiterScheme<Type,MeshType>> New
    (
        const fvMesh&,
        const word
    );

    // Selector from stream
    static autoPtr<limiterScheme<Type,MeshType>> New
    (
        const fvMesh&,
        Istream&
    );

    // Destructor
    virtual ~limiterScheme();

    // Compute the ratio of successive gradients at each face
    virtual tmp<meshField<faceScalar,MeshType>> r
    (
        const meshField<faceScalar,MeshType>& phi,
        const meshField<Type,MeshType>& field
    );

    // Return the limiter field
    virtual tmp<meshField<faceScalar,MeshType>> psi
    (
        const meshField<faceScalar,MeshType>& phi,
        const meshField<Type,MeshType>& field
    ) = 0;
};

}

}

}

#ifdef NoRepository
    #include "limiterScheme.C"
#endif

#endif
