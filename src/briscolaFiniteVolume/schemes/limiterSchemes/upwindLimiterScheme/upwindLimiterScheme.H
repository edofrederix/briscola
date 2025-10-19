#ifndef upwindLimiterScheme_H
#define upwindLimiterScheme_H

#include "limiterScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class upwindLimiterScheme
:
    public limiterScheme<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("linear");

    // Construct from mesh and stream
    upwindLimiterScheme(const fvMesh&, Istream& is);

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
    #include "upwindLimiterScheme.C"
#endif

#endif
