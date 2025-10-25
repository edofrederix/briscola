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
public:

    // Number of stencil components
    static const direction N = SType::nCsComponents;

    // Base smoother class

    class smoother
    {
    public:

        // Number of ghost cells
        static const label G = meshDirection<Type,MeshType>::ghosts;

        // Number of stencil components
        static const direction N = SType::nCsComponents;

        // Diagonal smooth at given level and direction
        static inline void smoothDiag
        (
            linearSystem<SType,Type,MeshType>& sys,
            const label l,
            const label d
        )
        {
            meshDirection<Type,MeshType>& x = sys.x()[l][d];

            const meshDirection<SType,MeshType>& A = sys.A()[l][d];
            const meshDirection<Type,MeshType>& b = sys.b()[l][d];
            const meshDirection<label,MeshType>& f = sys.forcingMask()[l][d];

            forAllCells(x, i, j, k)
                if (!f(i,j,k))
                    x(i,j,k) = b(i,j,k)/A(i,j,k).center();
        }

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
            const label l,
            const label sweeps,
            const labelList& converged
        ) = 0;
    };

    // Base external solver class

    class externalSolver
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
        TypeName("externalSolver");

        // Declare run-time constructor selection tables
        declareRunTimeSelectionTable
        (
            autoPtr,
            externalSolver,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& fvMsh,
                const label l
            ),
            (dict, fvMsh, l)
        );

        // Construct from dictionary, mesh and level
        externalSolver(const dictionary&, const fvMesh&, const label);

        // Selector from name and dictionary
        static autoPtr<externalSolver> New
        (
            const word&,
            const dictionary&,
            const fvMesh&,
            const label l = 0
        );

        // Selector from dictionary
        static autoPtr<externalSolver> New
        (
            const dictionary&,
            const fvMesh&,
            const label l = 0
        );

        // Destructor
        virtual ~externalSolver();

        // Access to the dictionary
        inline const dictionary& dict() const { return dict_; }

        // Access to the mesh
        inline const fvMesh& fvMsh() const { return fvMsh_; }

        // Return the mesh level number
        inline label levelNum() const { return l_; }

        // Return prepared switch
        inline Switch prepared() const { return prepared_; }

        // Prepare the solution of the given system
        virtual void prepare(linearSystem<SType,Type,MeshType>&)
        {
            prepared_ = true;
        }

        // Solve the given system
        virtual void solve(linearSystem<SType,Type,MeshType>&) = 0;
    };

protected:

    // Copy of the dictionary
    dictionary dict_;

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Smoother
    autoPtr<smoother> smoothPtr_;

    // Coarse grid solver
    autoPtr<externalSolver> coarseSolverPtr_;

    // Relative tolerance
    scalar relTol_;

    // Absolute tolerance
    scalar tolerance_;

    // Minimum number of iterations (defaults to 0)
    label minIter_;

    // Maximum number of iterations (defaults to 999)
    label maxIter_;

public:

    // Compute residual normalization factors of the system Ax = b at given mesh
    // level. The normalization factor is defined as sum(|Ax-Ay| + |b-Ay|), with
    // y = mean(x). The function reuses a pre-computed residual.
    static List<Type> normFactors
    (
        linearSystem<SType,Type,MeshType>&,
        const meshLevel<Type,MeshType>&
    );

    // For each direction, check if the solution is converged, given a current
    // and initial residual. Scalar residuals are computed as the
    // L-infinity-normal of Type
    static labelList checkConvergence
    (
        const List<Type>&,
        const List<Type>&,
        const scalar,
        const scalar
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

    // Solve the given system with optional constant matrix indicator
    virtual void solve
    (
        linearSystem<SType,Type,MeshType>&,
        const bool constMatrix = false
    ) = 0;

    // Solve the given tmp system
    virtual void solve
    (
        const tmp<linearSystem<SType,Type,MeshType>>&,
        const bool constMatrix = false
    );

    // Same as the static function but with default arguments
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

    // Set a singularity constraint at the given level
    void setSingularityConstraint
    (
        linearSystem<SType,Type,MeshType>& sys,
        const label l
    ) const;
};

}

}

}

#ifdef NoRepository
    #include "solver.C"
    #include "solverSmoother.C"
    #include "solverExternalSolver.C"
#endif

#endif
