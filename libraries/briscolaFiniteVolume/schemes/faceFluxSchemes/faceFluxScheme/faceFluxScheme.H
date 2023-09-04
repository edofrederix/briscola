#ifndef faceFluxScheme_H
#define faceFluxScheme_H

#include "scheme.H"
#include "FaceSpace.H"
#include "meshFields.H"
#include "meshLevels.H"
#include "meshDirections.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class faceFluxScheme
:
    public scheme
{
public:

    // Runtime type information
    TypeName("faceFluxScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        faceFluxScheme,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& fvMsh
        ),
        (dict, fvMsh)
    );

    // Construct from dictionary and mesh
    faceFluxScheme(const dictionary&, const fvMesh&);

    // Copy constructor
    faceFluxScheme(const faceFluxScheme&);

    // Selector
    static autoPtr<faceFluxScheme> New(const word, const fvMesh&);

    // Destructor
    virtual ~faceFluxScheme();

    // FaceFluxScheme of a colocated vector field
    virtual tmp<meshField<faceScalar,colocated>> faceFlux
    (
        const meshField<vector,colocated>&
    ) = 0;

    // FaceFluxScheme of a staggered scalar field
    virtual tmp<meshField<faceScalar,staggered>> faceFlux
    (
        const meshField<scalar,staggered>&
    ) = 0;
};

}

}

}

#endif
