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
    // The IB mask is equal to 0 in cell centers outside the IB,
    // and equal to 1 in cell centers inside the IB.
    meshField<label,MeshType> mask_;

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

} // end namespace fv

} // end namespace briscola

} // end namespace Foam

#ifdef NoRepository
    #include "penalization.C"
#endif

#endif
