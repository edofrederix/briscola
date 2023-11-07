#ifndef Mittal_H
#define Mittal_H

#include "immersedBoundaryMethod.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class Mittal
:
    public immersedBoundaryMethod<Type,MeshType>
{

public:

    // Runtime type information
    TypeName("Mittal");

    // Construct
    Mittal
    (
        dictionary&,
        const fvMesh&
    );

    // Destructor
    virtual ~Mittal();

    // Correct Jacobi points
    void correctJacobiPoints(meshLevel<Type,MeshType>&);

};

} // end namespace fv

} // end namespace briscola

} // end namespace Foam

#ifdef NoRepository
    #include "Mittal.C"
#endif

#endif
