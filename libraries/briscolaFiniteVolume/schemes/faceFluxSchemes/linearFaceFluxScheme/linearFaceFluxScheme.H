#ifndef linearFaceFluxScheme_H
#define linearFaceFluxScheme_H

#include "faceFluxScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class linearFaceFluxScheme
:
    public faceFluxScheme
{
public:

    // Runtime type information
    TypeName("linear");

    // Construct from dictionary and mesh
    linearFaceFluxScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    linearFaceFluxScheme(const fvMesh&);

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
    #include "linearFaceFluxScheme.C"
#endif

#endif
