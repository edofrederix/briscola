#ifndef linearInterpolationScheme_H
#define linearInterpolationScheme_H

#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class linearInterpolationScheme
:
    public interpolationScheme<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("linear");

    // Construct from stream and mesh
    linearInterpolationScheme(const fvMesh&, Istream& = scheme::nullStream);

    // Interpolation to the face
    virtual tmp<faceField<Type,MeshType>> interp
    (
        const meshField<Type,MeshType>&
    );
};

}

}

}

#ifdef NoRepository
    #include "linearInterpolationScheme.C"
#endif

#endif
