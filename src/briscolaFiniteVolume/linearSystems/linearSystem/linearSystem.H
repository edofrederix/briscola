#ifndef linearSystem_H
#define linearSystem_H

#include "linearSystemFwd.H"
#include "linearSystemAggregation.H"

#include "fvMesh.H"
#include "meshFields.H"
#include "OFstream.H"
#include "cachedRefCount.H"
#include "restrict.H"

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
    public cachedRefCount
{
public:

    // Number of ghost cells
    static const label G = meshDirection<Type,MeshType>::ghosts;

    // Number of stencil components
    static const direction N = SType::nCsComponents;

private:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Pointer to the unknown
    meshField<Type,MeshType>* xPtr_;

    // The stencil matrix which is deep by default
    meshField<SType,MeshType> A_;

    // Right-hand side source which is shallow by default
    meshField<Type,MeshType> b_;

    // Symmetry, diagonal and eliminated bools
    bool symmetric_;
    bool diagonal_;
    bool eliminated_;

    // IBM forcing mask
    autoPtr<meshField<label,MeshType>> forcingMask_;

    // Transfer data from other system
    void transfer(linearSystem<SType,Type,MeshType>&);

    // Compute residual for a direction
    template<bool CorrBcs = true, bool Mask = true>
    void residual(meshDirection<Type,MeshType>&) const;

    // Compute evaluation for a direction
    template<bool CorrBcs = true, bool Mask = true>
    void evaluate(meshDirection<Type,MeshType>&) const;

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

    // Construct from IOobject and mesh, without unknown
    linearSystem
    (
        const IOobject& io,
        const fvMesh& fvMsh
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

    // Temporary linear system constructor
    static tmp<linearSystem<SType,Type,MeshType>> New
    (
        const word& name,
        meshField<Type,MeshType>& x,
        const bool cached = false
    )
    {
        tmp<linearSystem<SType,Type,MeshType>> tSys =
            x.fvMsh().template
            NewObject<linearSystem<SType,Type,MeshType>>(name, cached);

        tSys->setUnknown(x);

        return tSys;
    }

    // Destructor
    virtual ~linearSystem();

    // Return the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Set the unknown
    inline void setUnknown(meshField<Type,MeshType>& x)
    {
        xPtr_ = &x;
    }

    // Return the unknown (unchecked!)
    inline meshField<Type,MeshType>& x() { return *xPtr_; }
    inline const meshField<Type,MeshType>& x() const { return *xPtr_; }

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

    // Return symmetry, diagonal and eliminated bools
    bool symmetric() const { return symmetric_; }
    bool diagonal() const { return diagonal_; }
    bool eliminated() const { return eliminated_; }

    // Access to the symmetry, diagonal eliminated bools
    bool& symmetric() { return symmetric_; }
    bool& diagonal() { return diagonal_; }
    bool& eliminated() { return eliminated_; }

    // A linear system is said to be singular if the matrix is symmetric,
    // non-diagonal and if the boundary conditions of the unknown are singular
    bool singular() const
    {
        return symmetric_ && !diagonal_ && xPtr_->operator[](0).singular();
    }

    // Compute the residual b-Ax and store in given object. Optionally don't
    // correct boundary conditions or apply no forcing mask. This is determined
    // via a template parameter instead of an argument, to avoid ambiguity with
    // other arguments.

    template<bool CorrBcs = true, bool Mask = true>
    void residual(meshField<Type,MeshType>&) const;

    template<bool CorrBcs = true, bool Mask = true>
    void residual(meshLevel<Type,MeshType>&) const;

    // Compute and return the residual b-Ax.Optionally don't correct boundary
    // conditions or apply no forcing mask. This is determined via a template
    // parameter instead of an argument, to avoid ambiguity with other
    // arguments.

    template<bool CorrBcs = true, bool Mask = true>
    tmp<meshField<Type,MeshType>> residual() const;

    template<bool CorrBcs = true, bool Mask = true>
    tmp<meshLevel<Type,MeshType>> residual(const label) const;

    // Compute the evaluation (Ax-b)/V and store in given object. Optionally
    // don't correct boundary conditions or apply no forcing mask. This is
    // determined via a template parameter instead of an argument, to avoid
    // ambiguity with other arguments.

    template<bool CorrBcs = true, bool Mask = true>
    void evaluate(meshField<Type,MeshType>&) const;

    template<bool CorrBcs = true, bool Mask = true>
    void evaluate(meshLevel<Type,MeshType>&) const;

    // Compute the evaluation (Ax-b)/V. Optionally don't correct boundary
    // conditions or apply no forcing mask. This is determined via a template
    // parameter instead of an argument, to avoid ambiguity with other
    // arguments.

    template<bool CorrBcs = true, bool Mask = true>
    tmp<meshField<Type,MeshType>> evaluate() const;

    template<bool CorrBcs = true, bool Mask = true>
    tmp<meshLevel<Type,MeshType>> evaluate(const label) const;

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
    void operator=(const zero);
    void operator=(const meshField<Type,MeshType>&);
    void operator=(const tmp<meshField<Type,MeshType>>&);

    void operator+=(const Type&);
    void operator+=(const List<Type>&);
    void operator+=(const meshField<Type,MeshType>&);
    void operator+=(const tmp<meshField<Type,MeshType>>&);
    template<class Type2> void operator+=(const Type2&);

    void operator-=(const Type&);
    void operator-=(const List<Type>&);
    void operator-=(const meshField<Type,MeshType>&);
    void operator-=(const tmp<meshField<Type,MeshType>>&);
    template<class Type2> void operator-=(const Type2&);

    // Mixed stencil type operations

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

    void operator*=(const scalar);
    void operator*=(const scalarList&);
    void operator*=(const meshField<scalar,MeshType>&);
    void operator*=(const tmp<meshField<scalar,MeshType>>&);

    void operator/=(const scalar);
    void operator/=(const scalarList&);
    void operator/=(const meshField<scalar,MeshType>&);
    void operator/=(const tmp<meshField<scalar,MeshType>>&);

    // Befriend solver
    template<class, class, class>
    friend class solver;
};

}

}

}

#include "linearSystemI.H"

#ifdef NoRepository
    #include "linearSystem.C"
#endif

#endif
