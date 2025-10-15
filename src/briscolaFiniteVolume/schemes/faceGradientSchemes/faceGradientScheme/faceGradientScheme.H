#ifndef faceGradientScheme_H
#define faceGradientScheme_H

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
        Istream,
        (
            const fvMesh& fvMsh,
            Istream& is
        ),
        (fvMsh, is)
    );

    // Construct from stream and mesh
    faceGradientScheme(const fvMesh&, Istream& = scheme::nullStream);

    // Copy constructor
    faceGradientScheme(const faceGradientScheme<Type,MeshType>&);

    // Selector
    static autoPtr<faceGradientScheme<Type,MeshType>> New
    (
        const fvMesh&,
        const word
    );

    // Destructor
    virtual ~faceGradientScheme();

    // Face-normal gradient
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
