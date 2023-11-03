#ifndef upwindFaceFluxScheme_H
#define upwindFaceFluxScheme_H

#include "faceFluxScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class upwindFaceFluxScheme
:
    public faceFluxScheme
{
public:

    // Runtime type information
    TypeName("upwind");

    // Construct from dictionary and mesh
    upwindFaceFluxScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    upwindFaceFluxScheme(const fvMesh&);

    // Face flux of a colocated scalar field
    virtual tmp<meshField<faceScalar,colocated>> faceFlux
    (
        const meshField<vector,colocated>&
    );

    // Face flux of a staggered scalar field
    virtual tmp<meshField<faceScalar,staggered>> faceFlux
    (
        const meshField<scalar,staggered>&
    );
};

}

}

}

#ifdef NoRepository
    #include "upwindFaceFluxScheme.C"
#endif

#endif
