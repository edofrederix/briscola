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

    // Construct from dictionary and mesh
    midPointInterpolationScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    midPointInterpolationScheme(const fvMesh&);

    // Interpolation of u to the face
    virtual tmp<meshField<Hex<Type>,MeshType>> interp
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
