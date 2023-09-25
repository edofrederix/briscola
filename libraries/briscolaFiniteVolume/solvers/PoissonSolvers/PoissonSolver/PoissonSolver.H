#ifndef PoissonSolver_H
#define PoissonSolver_H

#include "PoissonSolverFwd.H"

#include "linearSystem.H"
#include "NamedEnum.H"
#include "colocatedFields.H"
#include "staggeredFields.H"
#include "dictionary.H"
#include "immersedBoundary.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class PoissonSolver
{
protected:

    // Copy of the scheme dictionary
    dictionary dict_;

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Pointer to immersed boundary
    immersedBoundary<Type,MeshType>* IB_;

public:

    // Runtime type information
    TypeName("PoissonSolver");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        PoissonSolver,
        dictionary,
        (
            const word PoissonSolverName,
            const dictionary& dict,
            const fvMesh& fvMsh,
            immersedBoundary<Type,MeshType>* IB
        ),
        (PoissonSolverName, dict, fvMsh, IB)
    );

    // Construct from dictionary and mesh
    PoissonSolver(const dictionary&, const fvMesh&, immersedBoundary<Type,MeshType>*);

    // Selector
    static autoPtr<PoissonSolver<SType,Type,MeshType>> New
    (
        const word,
        const fvMesh&,
        immersedBoundary<Type,MeshType>*
    );

    // Destructor
    virtual ~PoissonSolver();

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Solve for the given unknown x, right-hand side source b and scalar
    // coefficient lambda the equation dx/dt = div(lambda*grad(x)) + b. Provide
    // lambda and b as pointers, so that they can be left empty which implies b
    // = 0 and lambda = 1. The final argument specifies whether or not to
    // include the d/dt term.

    virtual void solve
    (
        meshField<Type,MeshType>&,
        const meshField<Type,MeshType>*,
        const meshField<faceScalar,MeshType>*,
        const bool = false
    ) = 0;

    // Single-argument support function

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const bool ddt = false
    )
    {
        this->solve(x, nullptr, nullptr, ddt);
    }

    // Two-argument support function

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const meshField<Type,MeshType>* bPtr,
        const bool ddt = false
    )
    {
        this->solve(x, bPtr, nullptr, ddt);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const meshField<Type,MeshType>& b,
        const bool ddt = false
    )
    {
        this->solve(x, &b, nullptr, ddt);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const tmp<meshField<Type,MeshType>>& tb,
        const bool ddt = false
    )
    {
        this->solve(x, &tb(), nullptr, ddt);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const meshField<faceScalar,MeshType>* lambdaPtr,
        const bool ddt = false
    )
    {
        this->solve(x, nullptr, lambdaPtr, ddt);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const meshField<faceScalar,MeshType>& lambda,
        const bool ddt = false
    )
    {
        this->solve(x, nullptr, &lambda, ddt);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const tmp<meshField<faceScalar,MeshType>>& tlambda,
        const bool ddt = false
    )
    {
        this->solve(x, nullptr, &tlambda(), ddt);
    }

    // Three-argument support function

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const meshField<Type,MeshType>& b,
        const meshField<faceScalar,MeshType>* lambdaPtr,
        const bool ddt = false
    )
    {
        this->solve(x, &b, lambdaPtr, ddt);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const tmp<meshField<Type,MeshType>>& tb,
        const meshField<faceScalar,MeshType>* lambdaPtr,
        const bool ddt = false
    )
    {
        this->solve(x, &tb(), lambdaPtr, ddt);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const meshField<Type,MeshType>* bPtr,
        const meshField<faceScalar,MeshType>& lambda,
        const bool ddt = false
    )
    {
        this->solve(x, bPtr, &lambda, ddt);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const meshField<Type,MeshType>& b,
        const meshField<faceScalar,MeshType>& lambda,
        const bool ddt = false
    )
    {
        this->solve(x, &b, &lambda, ddt);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const tmp<meshField<Type,MeshType>>& tb,
        const meshField<faceScalar,MeshType>& lambda,
        const bool ddt = false
    )
    {
        this->solve(x, &tb(), &lambda, ddt);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const meshField<Type,MeshType>* bPtr,
        const tmp<meshField<faceScalar,MeshType>>& tlambda,
        const bool ddt = false
    )
    {
        this->solve(x, bPtr, &tlambda(), ddt);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const meshField<Type,MeshType>& b,
        const tmp<meshField<faceScalar,MeshType>>& tlambda,
        const bool ddt = false
    )
    {
        this->solve(x, &b, &tlambda(), ddt);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const tmp<meshField<Type,MeshType>>& tb,
        const tmp<meshField<faceScalar,MeshType>>& tlambda,
        const bool ddt = false
    )
    {
        this->solve(x, &tb(), &tlambda(), ddt);
    }
};

}

}

}

#ifdef NoRepository
    #include "PoissonSolver.C"
#endif

#endif
