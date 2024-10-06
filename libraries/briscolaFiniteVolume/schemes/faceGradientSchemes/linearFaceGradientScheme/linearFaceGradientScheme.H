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

    // Construct from stream and mesh
    linearFaceGradientScheme(const fvMesh&, Istream& = scheme::nullStream);

    // Face-normal gradient
    virtual tmp<meshField<LowerFaceSpace<Type>,MeshType>> faceGrad
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
