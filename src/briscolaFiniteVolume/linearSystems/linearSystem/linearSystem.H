#ifndef linearSystem_H
#define linearSystem_H

#include "linearSystemFwd.H"
#include "linearSystemAggregation.H"

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
    public regIOobject,
    public refCount
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

    // IBM forcing mask
    autoPtr<meshField<label,MeshType>> forcingMask_;

    // Transfer data from other system
    void transfer(linearSystem<SType,Type,MeshType>&);

public:

    // Runtime type information
    TypeName("linearSystem");

    // Construct from unknown
    linearSystem
    (
        const word name,
        meshField<Type,MeshType>& x,
        const bool registerObject = false
    );

    // Copy constructors
    linearSystem
    (
        const linearSystem<SType,Type,MeshType>& tSys,
        const bool registerObject = false
    );

    // Copy constructor from tmp
    linearSystem
    (
        const tmp<linearSystem<SType,Type,MeshType>>& tSys,
        const bool registerObject = false
    );

    // Copy constructors with new unknown
    linearSystem
    (
        const linearSystem<SType,Type,MeshType>& sys,
        meshField<Type,MeshType>& x,
        const bool registerObject = false
    );

    // Copy constructors from tmp with new unknown
    linearSystem
    (
        const tmp<linearSystem<SType,Type,MeshType>>& tSys,
        meshField<Type,MeshType>& x,
        const bool registerObject = false
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
    inline meshField<label,MeshType>& forcingMask()
    {
        return forcingMask_();
    }

    // Return const access to the IBM forcing mask
    inline const meshField<label,MeshType>& forcingMask() const
    {
        return forcingMask_();
    }

    // A linear system is said to be singular if all matrix row sums at the
    // coarsest level  are zero. See Remark 5.6.5 in "Multigrid", U. Trottenberg
    // et al., (2000). Result is cached.
    List<bool> singular(const bool clearCache = false);

    // Check if the linear system is diagonal. Result is cached.
    List<bool> diagonal(const bool clearCache = false);

    // Compute the residual b-Ax and store in given object. Optionally apply no
    // forcing mask. This is determined via a template parameter instead of an
    // argument, to avoid ambiguity with other arguments.

    template<bool NoMask = false>
    void residual(meshField<Type,MeshType>&) const;

    template<bool NoMask = false>
    void residual(meshLevel<Type,MeshType>&) const;

    template<bool NoMask = false>
    void residual(meshDirection<Type,MeshType>&) const;

    // Compute and return the residual b-Ax. Optionally apply no forcing mask.
    // This is determined via a template parameter instead of an argument, to
    // avoid ambiguity with other arguments.

    template<bool NoMask = false>
    tmp<meshField<Type,MeshType>> residual() const;

    template<bool NoMask = false>
    tmp<meshLevel<Type,MeshType>> residual(const label) const;

    template<bool NoMask = false>
    tmp<meshDirection<Type,MeshType>> residual(const label, const label) const;

    // Compute the evaluation (Ax-b)/V and store in given object. Optionally
    // apply no forcing mask. This is determined via a template parameter
    // instead of an argument, to avoid ambiguity with other arguments.

    template<bool NoMask = false>
    void evaluate(meshField<Type,MeshType>&) const;

    template<bool NoMask = false>
    void evaluate(meshLevel<Type,MeshType>&) const;

    template<bool NoMask = false>
    void evaluate(meshDirection<Type,MeshType>&) const;

    // Compute the evaluation (Ax-b)/V. Optionally apply no forcing mask. This
    // is determined via a template parameter instead of an argument, to avoid
    // ambiguity with other arguments.

    template<bool NoMask = false>
    tmp<meshField<Type,MeshType>> evaluate() const;

    template<bool NoMask = false>
    tmp<meshLevel<Type,MeshType>> evaluate(const label) const;

    template<bool NoMask = false>
    tmp<meshDirection<Type,MeshType>> evaluate(const label, const label) const;

    // Eliminate ghosts from the linear system where possible
    void eliminateGhosts();

    // Set the IBM forcing mask if not yet done
    void setForcingMask();

    // Allow the IBM forcing mask to be set on a const linear system
    void setForcingMask() const
    {
        const_cast<linearSystem<SType,Type,MeshType>*>(this)->setForcingMask();
    }

    // Write data (not supported)
    bool writeData(Ostream& os) const { return os; }

    // Write the system at given level to an output stream. The last column is
    // the right-hand side source. The dense matrix is written, so use only on
    // very coarse meshes for debugging.
    bool writeLevel(Ostream& os, const label = 0) const;

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

}

}

}

#include "linearSystemI.H"

#ifdef NoRepository
    #include "linearSystem.C"
#endif

#endif
