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

    // Construct from dictionary and mesh
    linearInterpolationScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    linearInterpolationScheme(const fvMesh&);

    // Interpolation to the face
    virtual tmp<meshField<LowerFaceSpace<Type>,MeshType>> interp
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
