#ifndef faceGradientScheme_H
#define faceGradientScheme_H

#include "scheme.H"
#include "FaceSpace.H"
#include "meshFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class faceGradientScheme
:
    public scheme
{
public:

    // Runtime type information
    TypeName("faceGradientScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        faceGradientScheme,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& fvMsh
        ),
        (dict, fvMsh)
    );

    // Construct from dictionary and mesh
    faceGradientScheme(const dictionary&, const fvMesh&);

    // Copy constructor
    faceGradientScheme(const faceGradientScheme<Type,MeshType>&);

    // Selector
    static autoPtr<faceGradientScheme<Type,MeshType>> New(const word, const fvMesh&);

    // Destructor
    virtual ~faceGradientScheme();

    // Face-normal gradient of u
    virtual tmp<meshField<FaceSpace<Type>,MeshType>> faceGrad
    (
        const meshField<Type,MeshType>&
    ) = 0;
};

}

}

}

#ifdef NoRepository
    #include "faceGradientScheme.C"
#endif

#endif
