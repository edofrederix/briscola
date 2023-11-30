#ifndef solver_H
#define solver_H

#include "solverFwd.H"

#include "linearSystem.H"
#include "colocatedFields.H"
#include "staggeredFields.H"
#include "dictionary.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class solver
{
protected:

    // Copy of the scheme dictionary
    dictionary dict_;

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Relative tolerance
    scalar relTol_;

    // Absolute tolerance
    scalar tolerance_;

    // Minimum number of iterations (defaults to 0)
    label minIter_;

    // Maximum number of iterations (defaults to 999)
    label maxIter_;

public:

    // Runtime type information
    TypeName("solver");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        solver,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& fvMsh
        ),
        (dict, fvMsh)
    );

    // Construct from dictionary and mesh
    solver(const dictionary&, const fvMesh&);

    // Selector from solver name
    static autoPtr<solver<SType,Type,MeshType>> New
    (
        const word,
        const fvMesh&
    );

    // Selector from dictionary
    static autoPtr<solver<SType,Type,MeshType>> New
    (
        const dictionary&,
        const fvMesh&
    );

    // Destructor
    virtual ~solver();

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Solve the given system
    virtual void solve(linearSystem<SType,Type,MeshType>&) = 0;

    // Compute residual normalization factors of the system Ax=b at given mesh
    // level. The normalization factor is defined as sum(|Ax-Ay|+|b-Ay|), with y
    // = mean(x)
    virtual List<Type> normFactors
    (
        const linearSystem<SType,Type,MeshType>&,
        const label
    ) const;

    // Red-black Gauss-Seidel smoother
    static void RBGS
    (
        linearSystem<SType,Type,MeshType>& sys,
        List<Type>& xi,
        const label l,
        const label sweeps,
        const labelList& converged,
        const scalar omega = 1.0
    );

    // Lexicographic Gauss-Seidel smoother
    static void LEXGS
    (
        linearSystem<SType,Type,MeshType>& sys,
        List<Type>& xi,
        const label l,
        const label sweeps,
        const labelList& converged,
        const scalar omega = 1.0
    );

    // Jacobi smoother
    static void JAC
    (
        linearSystem<SType,Type,MeshType>& sys,
        List<Type>& xi,
        const label l,
        const label sweeps,
        const labelList& converged,
        const scalar omega = 0.8
    );

    // Print the statistics of the solver outcome for given solver name, field
    // name, list of initial residual, final residuals and number of iterations
    static void printSolverStats
    (
        const word,
        const word,
        const List<Type>,
        const List<Type>,
        const label
    );

    // For each direction, check if the solution is converged, given a current
    // and initial residual. Scalar residuals are computed as the
    // L-infinity-normal of Type
    labelList checkConvergence
    (
        const List<Type>&,
        const List<Type>&,
        const scalar,
        const scalar
    ) const;

    // Same function with default arguments
    labelList checkConvergence
    (
        const List<Type>& currentResiduals,
        const List<Type>& initialResiduals
    ) const
    {
        return
            checkConvergence
            (
                currentResiduals,
                initialResiduals,
                this->relTol_,
                this->tolerance_
            );
    }
};

}

}

}

#ifdef NoRepository
    #include "solver.C"
#endif

#endif
