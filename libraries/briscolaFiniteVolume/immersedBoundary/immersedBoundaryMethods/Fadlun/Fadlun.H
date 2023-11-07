#ifndef Fadlun_H
#define Fadlun_H

#include "immersedBoundaryMethod.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class Fadlun
:
    public immersedBoundaryMethod<Type,MeshType>
{
    // The wall adjacent mask is equal to 1 in fluid cells with at
    // least one neighboring cell in the IB, and 0 in other cells
    meshField<label,MeshType> wallAdjMask_;

    // The wall distance field stores the wall distances to the
    // IB in the components of the stencils
    meshField<faceScalar,MeshType> wallDist_;

public:

    // Runtime type information
    TypeName("Fadlun");

    // Construct
    Fadlun
    (
        dictionary&,
        const fvMesh&
    );

    // Destructor
    virtual ~Fadlun();

    // Correct linear system
    void correctLinearSystem(linearSystem<stencil,Type,MeshType>&);

};

}

}

}

#ifdef NoRepository
    #include "Fadlun.C"
#endif

#endif
