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

    // Singularity bools
    List<bool> singular_;

    // Diagonal bools
    List<bool> diagonal_;

    // IBM forcing points mask
    meshField<label,MeshType> IBMForcingMask_;

    // Can this linear system be modified by an IBM
    const bool IBM_;

    // Transfer data from other system
    void transfer(linearSystem<SType,Type,MeshType>&);

    // Set IBM forcing mask
    void setIBMForcingMask();

public:

    // Runtime type information
    TypeName("linearSystem");

    // Construct from unknown
    linearSystem(meshField<Type,MeshType>& x, bool IBM = false);

    // Copy constructors
    linearSystem(const linearSystem<SType,Type,MeshType>&);
    linearSystem(const tmp<linearSystem<SType,Type,MeshType>>&);

    // Copy constructors with new unknown
    linearSystem
    (
        const linearSystem<SType,Type,MeshType>& sys,
        meshField<Type,MeshType>& x,
        bool IBM = false
    );

    linearSystem
    (
        const tmp<linearSystem<SType,Type,MeshType>>& tSys,
        meshField<Type,MeshType>& x,
        bool IBM = false
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

    // Return the IBM forcing mask
    inline meshField<label,MeshType>& IBMForcingMask() { return IBMForcingMask_; }
    inline const meshField<label,MeshType>& IBMForcingMask() const { return IBMForcingMask_; }

    // A linear system is said to be singular if all matrix row sums at the
    // coarsest level  are zero. See Remark 5.6.5 in "Multigrid", U. Trottenberg
    // et al., (2000). Result is cached.
    List<bool> singular(const bool clearCache = false);

    // Check if the linear system is diagonal. Result is cached.
    List<bool> diagonal(const bool clearCache = false);

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
    void operator+=(const List<Type>&);
    void operator+=(const meshField<Type,MeshType>&);
    void operator+=(const tmp<meshField<Type,MeshType>>&);

    void operator-=(const Type&);
    void operator-=(const List<Type>&);
    void operator-=(const meshField<Type,MeshType>&);
    void operator-=(const tmp<meshField<Type,MeshType>>&);

    template<class SType2>
    void operator=(const linearSystem<SType2,Type,MeshType>&);
    template<class SType2>
    void operator=(const tmp<linearSystem<SType2,Type,MeshType>>&);

    template<class SType2>
    void operator+=(const SType2&);
    template<class SType2>
    void operator+=(const linearSystem<SType2,Type,MeshType>&);
    template<class SType2>
    void operator+=(const tmp<linearSystem<SType2,Type,MeshType>>&);

    template<class SType2>
    void operator-=(const SType2&);
    template<class SType2>
    void operator-=(const linearSystem<SType2,Type,MeshType>&);
    template<class SType2>
    void operator-=(const tmp<linearSystem<SType2,Type,MeshType>>&);

    void operator*=(const scalar);
    void operator*=(const scalarList&);
    void operator*=(const meshField<scalar,MeshType>&);
    void operator*=(const tmp<meshField<scalar,MeshType>>&);

    void operator/=(const scalar);
    void operator/=(const scalarList&);
    void operator/=(const meshField<scalar,MeshType>&);
    void operator/=(const tmp<meshField<scalar,MeshType>>&);
};

// Write a stencil system at given level to a file. The last column is the
// right-hand side source. Full matrix is written, so use only on very coarse
// meshes.
template<class SType, class Type, class MeshType>
void writeToFile
(
    const linearSystem<SType,Type,MeshType>&,
    const fileName file,
    const label = 0
);

}

}

}

#include "linearSystemI.H"

#ifdef NoRepository
    #include "linearSystem.C"
#endif

#endif
