#ifndef interpolationScheme_H
#define interpolationScheme_H

#include "scheme.H"
#include "meshFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class interpolationScheme
:
    public scheme
{
public:

    // Runtime type information
    TypeName("interpolationScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        interpolationScheme,
        Istream,
        (
            const fvMesh& fvMsh,
            Istream& is
        ),
        (fvMsh, is)
    );

    // Construct from stream and mesh
    interpolationScheme(const fvMesh&, Istream& = scheme::nullStream);

    // Copy constructor
    interpolationScheme(const interpolationScheme<Type,MeshType>&);

    // Selector
    static autoPtr<interpolationScheme<Type,MeshType>> New
    (
        const fvMesh&,
        const word
    );

    // Destructor
    virtual ~interpolationScheme();

    // Interpolation to the face
    virtual tmp<meshField<LowerFaceSpace<Type>,MeshType>> interp
    (
        const meshField<Type,MeshType>&
    ) = 0;

    virtual tmp<meshField<LowerFaceSpace<Type>,MeshType>> interp
    (
        const tmp<meshField<Type,MeshType>>& tField
    )
    {
        tmp<meshField<LowerFaceSpace<Type>,MeshType>> tInterp
        (
            interp(tField())
        );

        if (tField.isTmp())
            tField.clear();

        return tInterp;
    }
};

// Interpolation functions from colocated to staggered

template<class Type>
tmp<meshField<Type,staggered>> stagInterp
(
    const meshField<Type,colocated>&
);

template<class Type>
tmp<meshField<Type,staggered>> stagInterp
(
    const tmp<meshField<Type,colocated>>&
);

// Interpolation functions from colocated to staggered faces

template<class Type>
tmp<meshField<LowerFaceSpace<Type>,staggered>> stagFaceInterp
(
    const meshField<Type,colocated>&
);

template<class Type>
tmp<meshField<LowerFaceSpace<Type>,staggered>> stagFaceInterp
(
    const tmp<meshField<Type,colocated>>&
);

// Interpolation functions from staggered to colocated

template<class Type>
tmp<meshField<Type,colocated>> coloInterp
(
    const meshField<Type,staggered>&
);

template<class Type>
tmp<meshField<Type,colocated>> coloInterp
(
    const tmp<meshField<Type,staggered>>&
);

// Interpolation function from staggered to colocated faces

template<class Type>
tmp<meshField<LowerFaceSpace<Type>,colocated>> colorFaceInterp
(
    const meshField<Type,staggered>&
);

template<class Type>
tmp<meshField<LowerFaceSpace<Type>,colocated>> colorFaceInterp
(
    const tmp<meshField<Type,staggered>>&
);

}

}

}

#ifdef NoRepository
    #include "interpolationScheme.C"
#endif

#endif
