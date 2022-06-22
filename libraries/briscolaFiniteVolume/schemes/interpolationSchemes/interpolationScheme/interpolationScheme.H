#ifndef interpolationScheme_H
#define interpolationScheme_H

#include "scheme.H"
#include "meshFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class interpolationScheme
:
    public scheme
{
public:

    // Runtime type information
    TypeName("interpolationScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        interpolationScheme,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& fvMsh
        ),
        (dict, fvMsh)
    );

    // Construct from dictionary and mesh
    interpolationScheme(const dictionary&, const fvMesh&);

    // Copy constructor
    interpolationScheme(const interpolationScheme<Type,MeshType>&);

    // Selector
    static autoPtr<interpolationScheme<Type,MeshType>> New(const word, const fvMesh&);

    // Destructor
    virtual ~interpolationScheme();

    // Interpolation of u to the face
    virtual tmp<meshField<Hex<Type>,MeshType>> interp
    (
        const meshField<Type,MeshType>&
    ) = 0;
};

}

}

}

#ifdef NoRepository
    #include "interpolationScheme.C"
#endif

#endif
