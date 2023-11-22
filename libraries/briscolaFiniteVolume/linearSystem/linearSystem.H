#ifndef linearSystem_H
#define linearSystem_H

#include "linearSystemFwd.H"

#include "fvMesh.H"
#include "meshFields.H"
#include "OFstream.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class linearSystem
:
    public tmp<linearSystem<SType,Type,MeshType>>::refCount
{
    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Reference to the unknown
    meshField<Type,MeshType>& x_;

    // The stencil matrix which is deep by default
    meshField<SType,MeshType> A_;

    // Right-hand side source which is shallow by default
    meshField<Type,MeshType> b_;

    // Bool that keeps track of system singularity
    bool singular_;

    // Transfer data from other system
    void transfer(linearSystem<SType,Type,MeshType>&);

public:

    // Runtime type information
    TypeName("linearSystem");

    // Construct from unknown
    linearSystem(meshField<Type,MeshType>&);

    // Copy constructors
    linearSystem(const linearSystem<SType,Type,MeshType>&);
    linearSystem(const tmp<linearSystem<SType,Type,MeshType>>&);

    // Copy constructors with new unknown
    linearSystem
    (
        const linearSystem<SType,Type,MeshType>&,
        meshField<Type,MeshType>&
    );

    linearSystem
    (
        const tmp<linearSystem<SType,Type,MeshType>>&,
        meshField<Type,MeshType>&
    );

    // Destructor
    virtual ~linearSystem();

    // Return the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Return the unknown
    inline meshField<Type,MeshType>& x() { return x_; }
    inline const meshField<Type,MeshType>& x() const { return x_; }

    // Return the stencil matrix
    inline meshField<SType,MeshType>& A() { return A_; }
    inline const meshField<SType,MeshType>& A() const { return A_; }

    // Return the right-hand side source
    inline meshField<Type,MeshType>& b() { return b_; }
    inline const meshField<Type,MeshType>& b() const { return b_; }

    // Return the singularity bool
    inline bool& singular() { return singular_; }
    inline const bool& singular() const { return singular_; }

    // Compute the residual b-Ax and store in given object
    void residual(meshField<Type,MeshType>&) const;
    void residual(meshLevel<Type,MeshType>&) const;
    void residual(meshDirection<Type,MeshType>&) const;

    // Compute the residual b-Ax
    tmp<meshField<Type,MeshType>> residual() const;
    tmp<meshLevel<Type,MeshType>> residual(const label) const;
    tmp<meshDirection<Type,MeshType>> residual(const label, const label) const;

    // Compute the evaluation (Ax-b)/V and store in given object
    void evaluate(meshField<Type,MeshType>&) const;
    void evaluate(meshLevel<Type,MeshType>&) const;
    void evaluate(meshDirection<Type,MeshType>&) const;

    // Compute the evaluation (Ax-b)/V
    tmp<meshField<Type,MeshType>> evaluate() const;
    tmp<meshLevel<Type,MeshType>> evaluate(const label) const;
    tmp<meshDirection<Type,MeshType>> evaluate(const label, const label) const;

    // Eliminate ghosts from the linear system where possible
    void eliminateGhosts();

    // Eliminate ghosts from the linear system where possible, at given level
    void eliminateGhosts(const label);

    // Member Operators

    void operator=(const linearSystem<SType,Type,MeshType>&);
    void operator=(const tmp<linearSystem<SType,Type,MeshType>>&);

    void operator+=(const linearSystem<SType,Type,MeshType>&);
    void operator+=(const tmp<linearSystem<SType,Type,MeshType>>&);

    void operator-=(const linearSystem<SType,Type,MeshType>&);
    void operator-=(const tmp<linearSystem<SType,Type,MeshType>>&);

    void operator=(const Type&);
    void operator=(const meshField<Type,MeshType>&);
    void operator=(const tmp<meshField<Type,MeshType>>&);

    void operator+=(const Type&);
    void operator+=(const meshField<Type,MeshType>&);
    void operator+=(const tmp<meshField<Type,MeshType>>&);

    void operator-=(const Type&);
    void operator-=(const meshField<Type,MeshType>&);
    void operator-=(const tmp<meshField<Type,MeshType>>&);

    template<class SType2>
    void operator=(const linearSystem<SType2,Type,MeshType>&);
    template<class SType2>
    void operator=(const tmp<linearSystem<SType2,Type,MeshType>>&);

    template<class SType2>
    void operator+=(const linearSystem<SType2,Type,MeshType>&);
    template<class SType2>
    void operator+=(const tmp<linearSystem<SType2,Type,MeshType>>&);

    template<class SType2>
    void operator-=(const linearSystem<SType2,Type,MeshType>&);
    template<class SType2>
    void operator-=(const tmp<linearSystem<SType2,Type,MeshType>>&);
};

// Write a stencil system at given level to a file. The last column is the
// right-hand side source. Full matrix is written, so use only on very coarse
// meshes.
template<class SType, class Type, class MeshType>
void writeToFile
(
    const linearSystem<SType,Type,MeshType>&,
    OFstream&,
    const label
);

}

}

}

#include "linearSystemI.H"

#ifdef NoRepository
    #include "linearSystem.C"
#endif

#endif
