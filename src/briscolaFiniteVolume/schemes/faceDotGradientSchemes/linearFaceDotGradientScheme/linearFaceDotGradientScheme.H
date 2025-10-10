#ifndef linearFaceDotGradientScheme_H
#define linearFaceDotGradientScheme_H

#include "faceDotGradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Forward declaration

template<class Type, class MeshType>
class linearFaceDotGradientScheme;

// Colocated specialization

template<class Type>
class linearFaceDotGradientScheme<Type,colocated>
:
    public faceDotGradientScheme<Type,colocated>
{
protected:

    // Generate and store cache fields
    virtual void cache(const fvMesh&) const;

public:

    // Runtime type information
    TypeName("linear");

    // Construct from stream and mesh
    linearFaceDotGradientScheme
    (
        const fvMesh& fvMsh,
        Istream& is = scheme::nullStream
    )
    :
        faceDotGradientScheme<Type,colocated>(fvMsh, is)
    {}

    // Face-dot gradient multiplied by face area
    virtual tmp<meshField<FaceSpace<Type>,colocated>> faceDotGrad
    (
        const meshField<Type,colocated>&
    );
};

// Staggered specialization

template<class Type>
class linearFaceDotGradientScheme<Type,staggered>
:
    public faceDotGradientScheme<Type,staggered>
{
public:

    // Runtime type information
    TypeName("linear");

    // Construct from stream and mesh
    linearFaceDotGradientScheme
    (
        const fvMesh& fvMsh,
        Istream& is = scheme::nullStream
    )
    :
        faceDotGradientScheme<Type,staggered>(fvMsh, is)
    {}

    // Face-normal inner product gradient
    virtual tmp<meshField<FaceSpace<Type>,staggered>> faceDotGrad
    (
        const meshField<Type,staggered>&
    );
};

}

}

}

#ifdef NoRepository
    #include "linearFaceDotGradientScheme.C"
#endif

#endif
