#ifndef solver_H
#define solver_H

#include "solverFwd.H"

#include "linearSystem.H"
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
class solver
{
public:

    // Base smoother class with singular system augmentation, see Sec. 5.6.4 of
    // "Multigrid", U. Trottenberg et al., 2000. If the xi list is empty, no
    // singular system augmentation is applied. Having xi as an argument enables
    // its reuse between different calls to the smoothers.

    class smoother
    {
    protected:

        // Reference to the solver dictionary
        const dictionary& dict_;

        // Reference to the mesh
        const fvMesh& fvMsh_;

    public:

        // Runtime type information
        TypeName("smoother");

        // Declare run-time constructor selection tables
        declareRunTimeSelectionTable
        (
            autoPtr,
            smoother,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& fvMsh
            ),
            (dict, fvMsh)
        );

        // Construct from dictionary and mesh
        smoother(const dictionary&, const fvMesh&);

        // Selector from name
        static autoPtr<smoother> New
        (
            const word&,
            const dictionary&,
            const fvMesh&
        );

        // Destructor
        virtual ~smoother();

        // Access to the dictionary
        inline const dictionary& dict() const { return dict_; }

        // Access to the mesh
        inline const fvMesh& fvMsh() const { return fvMsh_; }

        // Smooth at given level
        virtual void smooth
        (
            linearSystem<SType,Type,MeshType>& sys,
            List<Type>& xi,
            const label l,
            const label sweeps,
            const labelList& converged
        ) = 0;
    };

    // Base direct solver class

    class directSolver
    {
    protected:

        // Reference to the solver dictionary
        const dictionary& dict_;

        // Reference to the mesh
        const fvMesh& fvMsh_;

        // Level number
        const label l_;

        // Prepared switch
        Switch prepared_;

    public:

        // Runtime type information
        TypeName("directSolver");

        // Declare run-time constructor selection tables
        declareRunTimeSelectionTable
        (
            autoPtr,
            directSolver,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& fvMsh,
                const label l
            ),
            (dict, fvMsh, l)
        );

        // Construct from dictionary and mesh
        directSolver(const dictionary&, const fvMesh&, const label);

        // Selector from name
        static autoPtr<directSolver> New
        (
            const word&,
            const dictionary&,
            const fvMesh&,
            const label l = 0
        );

        // Destructor
        virtual ~directSolver();

        // Access to the dictionary
        inline const dictionary& dict() const { return dict_; }

        // Access to the mesh
        inline const fvMesh& fvMsh() const { return fvMsh_; }

        // Return the mesh level number
        inline label levelNum() const { return l_; }

        // Return prepared switch
        inline Switch prepared() const { return prepared_; }

        // Prepare the solution of the given system, indicating if singular
        // system augmentation should be applied
        virtual void prepare
        (
            const linearSystem<SType,Type,MeshType>&,
            const List<bool>&
        )
        {
            prepared_ = true;
        }

        // Solve the given system, indicating if singular system augmentation
        // should be applied
        virtual void solve
        (
            linearSystem<SType,Type,MeshType>&,
            const List<bool>&
        ) = 0;
    };

protected:

    // Copy of the dictionary
    dictionary dict_;

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Smoother
    autoPtr<smoother> smoothPtr_;

    // Direct solver
    autoPtr<directSolver> directSolvePtr_;

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
    solver
    (
        const dictionary&,
        const fvMesh&
    );

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

    // Solve the given system with optional constant matrix indicator
    virtual void solve
    (
        linearSystem<SType,Type,MeshType>&,
        const bool constMatrix = false
    ) = 0;

    // Compute residual normalization factors of the system Ax=b at given mesh
    // level. The normalization factor is defined as sum(|Ax-Ay|+|b-Ay|), with y
    // = mean(x). The function reuses the residual that should be pre-computed,
    // so that it can be re-used.
    virtual List<Type> normFactors
    (
        const linearSystem<SType,Type,MeshType>&,
        const meshField<Type,MeshType>&,
        const label
    ) const;

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
    #include "solverSmoother.C"
    #include "solverDirectSolver.C"
#endif

#endif
