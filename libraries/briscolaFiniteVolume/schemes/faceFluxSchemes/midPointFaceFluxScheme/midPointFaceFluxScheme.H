#ifndef midPointFaceFluxScheme_H
#define midPointFaceFluxScheme_H

#include "faceFluxScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class midPointFaceFluxScheme
:
    public faceFluxScheme
{
public:

    // Runtime type information
    TypeName("midPoint");

    // Construct from dictionary and mesh
    midPointFaceFluxScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    midPointFaceFluxScheme(const fvMesh&);

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
    #include "midPointFaceFluxScheme.C"
#endif

#endif
