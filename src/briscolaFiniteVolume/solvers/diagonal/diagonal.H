#ifndef diagonal_H
#define diagonal_H

#include "solver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class diagonal
:
    public solver<SType,Type,MeshType>
{
public:

    // Runtime type information
    TypeName("diagonal");

    // Construct from dictionary and mesh
    diagonal(const dictionary&, const fvMesh&);

    // Solve the given system with optional constant matrix indicator
    virtual void solve
    (
        linearSystem<SType,Type,MeshType>& sys,
        const bool constMatrix = false
    );
};

}

}

}

#ifdef NoRepository
    #include "diagonal.C"
#endif

#endif
