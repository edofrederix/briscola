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
    tmp<faceField<scalar,MeshType>> rType
    (
        const faceField<scalar,MeshType>& phi,
        const meshField<scalar,MeshType>& field,
        const bool deep
    );

    // Vector helper function for the ratio of successive gradients
    tmp<faceField<scalar,MeshType>> rType
    (
        const faceField<scalar,MeshType>& phi,
        const meshField<vector,MeshType>& field,
        const bool deep
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
    virtual tmp<faceField<scalar,MeshType>> r
    (
        const faceField<scalar,MeshType>& phi,
        const meshField<Type,MeshType>& field,
        const bool deep
    );

    // Return the limiter field
    virtual tmp<faceField<scalar,MeshType>> psi
    (
        const faceField<scalar,MeshType>& phi,
        const meshField<Type,MeshType>& field,
        const bool deep
    ) = 0;
};

}

}

}

#ifdef NoRepository
    #include "limiterScheme.C"
#endif

#endif
