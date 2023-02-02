#ifndef FFTPoissonSolver_H
#define FFTPoissonSolver_H

#include "PoissonSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class FFTPoissonSolver
:
    public PoissonSolver<SType,Type,MeshType>
{
public:

    // Runtime type information
    TypeName("FFT");

    // Construct from dictionary and mesh
    FFTPoissonSolver
    (
        const word,
        const dictionary&,
        const fvMesh&
    );

    // Solve function
    virtual void solve
    (
        meshField<Type,MeshType>&,
        const meshField<Type,MeshType>*,
        const meshField<faceScalar,MeshType>*,
        const bool = false
    );
};

}

}

}

#ifdef NoRepository
    #include "FFTPoissonSolver.C"
#endif

#endif
