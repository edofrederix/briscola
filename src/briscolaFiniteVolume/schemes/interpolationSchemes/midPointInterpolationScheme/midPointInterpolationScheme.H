#ifndef midPointInterpolationScheme_H
#define midPointInterpolationScheme_H

#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class midPointInterpolationScheme
:
    public interpolationScheme<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("midPoint");

    // Construct from stream and mesh
    midPointInterpolationScheme(const fvMesh&, Istream& = scheme::nullStream);

    // Interpolation to the face
    virtual tmp<meshField<FaceSpace<Type>,MeshType>> interp
    (
        const meshField<Type,MeshType>&
    );
};

}

}

}

#ifdef NoRepository
    #include "midPointInterpolationScheme.C"
#endif

#endif
