#ifndef midPointFaceGradientScheme_H
#define midPointFaceGradientScheme_H

#include "faceGradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class midPointFaceGradientScheme
:
    public faceGradientScheme<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("midPoint");

    // Construct from dictionary and mesh
    midPointFaceGradientScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    midPointFaceGradientScheme(const fvMesh&);

    // Face-normal gradient of u
    virtual tmp<meshField<Hex<Type>,MeshType>> faceGrad
    (
        const meshField<Type,MeshType>&
    );
};

}

}

}

#ifdef NoRepository
    #include "midPointFaceGradientScheme.C"
#endif

#endif
