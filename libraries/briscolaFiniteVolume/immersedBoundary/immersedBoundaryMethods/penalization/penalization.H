#ifndef penalization_H
#define penalization_H

#include "immersedBoundaryMethod.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class penalization
:
    public immersedBoundaryMethod<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("penalization");

    // Construct
    penalization
    (
        dictionary&,
        const fvMesh&
    );

    // Destructor
    virtual ~penalization();

    // Correct linear system
    void correctLinearSystem(linearSystem<stencil,Type,MeshType>&);

};

}

}

}

#ifdef NoRepository
    #include "penalization.C"
#endif

#endif
