#ifndef curvatureInterpolationScheme_H
#define curvatureInterpolationScheme_H

#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class curvatureInterpolationScheme
:
    public interpolationScheme<Type,MeshType>
{

    scalar threshold_;

public:

    // Runtime type information
    TypeName("curvature");

    // Construct from dictionary and mesh
    curvatureInterpolationScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    curvatureInterpolationScheme(const fvMesh&);

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
    #include "curvatureInterpolationScheme.C"
#endif

#endif
