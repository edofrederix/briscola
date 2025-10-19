#ifndef twoPhaseLimiterScheme_H
#define twoPhaseLimiterScheme_H

#include "limiterScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class twoPhaseLimiterScheme
:
    public limiterScheme<Type,MeshType>
{
    // The actual limiter scheme
    autoPtr<limiterScheme<Type,MeshType>> limiter_;

public:

    // Runtime type information
    TypeName("twoPhase");

    // Construct from mesh and stream
    twoPhaseLimiterScheme(const fvMesh&, Istream& is);

    // Return the limiter field
    virtual tmp<faceField<scalar,MeshType>> psi
    (
        const faceField<scalar,MeshType>& phi,
        const meshField<Type,MeshType>& field
    );
};

}

}

}

#ifdef NoRepository
    #include "twoPhaseLimiterScheme.C"
#endif

#endif
