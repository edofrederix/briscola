#ifndef linearFaceGradientScheme_H
#define linearFaceGradientScheme_H

#include "faceGradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class linearFaceGradientScheme
:
    public faceGradientScheme<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("linear");

    // Construct from dictionary and mesh
    linearFaceGradientScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    linearFaceGradientScheme(const fvMesh&);

    // Face-normal gradient of u
    virtual tmp<meshField<FaceSpace<Type>,MeshType>> faceGrad
    (
        const meshField<Type,MeshType>&
    );
};

}

}

}

#ifdef NoRepository
    #include "linearFaceGradientScheme.C"
#endif

#endif
