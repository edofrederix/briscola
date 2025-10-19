#ifndef PoissonSolver_H
#define PoissonSolver_H

#include "PoissonSolverFwd.H"

#include "linearSystem.H"
#include "NamedEnum.H"
#include "colocatedFields.H"
#include "staggeredFields.H"
#include "dictionary.H"
#include "RungeKuttaScheme.H"
#include "FastPtrList.H"

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

    // Flux field
    autoPtr<faceField<Type,MeshType>> fluxPtr_;

    // Switch to compute the flux (defaults to false)
    Switch computeFlux_;

    // Pointer to the Runge-Kutta scheme, if found
    RungeKuttaScheme* rkSchemePtr_;

    // Runge-Kutta stage solution cache which can be used by derived types as
    // initial guess to improve convergence
    FastPtrList<meshField<Type,MeshType>> stages_;

    // Function to initialize the flux field
    virtual void initFlux();

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
            const dictionary& dict,
            const fvMesh& fvMsh
        ),
        (dict, fvMsh)
    );

    // Construct from dictionary and mesh
    PoissonSolver(const dictionary&, const fvMesh&);

    // Selector from dictionary
    static autoPtr<PoissonSolver<SType,Type,MeshType>> New
    (
        const dictionary&,
        const fvMesh&
    );

    // Selector from solver name
    static autoPtr<PoissonSolver<SType,Type,MeshType>> New
    (
        const word,
        const fvMesh&
    );

    // Destructor
    virtual ~PoissonSolver();

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Return the Runge-Kutta scheme pointer
    inline RungeKuttaScheme* rkSchemePtr() const
    {
        return rkSchemePtr_;
    }

    // Access to the Runge-Kutta stage cache
    inline FastPtrList<meshField<Type,MeshType>>& stages()
    {
        return stages_;
    }

    // Return the Runge-Kutta stage cache
    inline const FastPtrList<meshField<Type,MeshType>>& stages() const
    {
        return stages_;
    }

    // Access to the flux
    inline const faceField<Type,MeshType>& flux() const
    {
        #ifdef FULLDEBUG

        if (fluxPtr_.empty())
        {
            FatalErrorInFunction
                << "Flux not set" << endl << abort(FatalError);
        }

        #endif

        return fluxPtr_();
    }

    // Computation of the flux switch
    inline Switch computeFlux() const
    {
        return computeFlux_;
    }

    // Enable the computation of the flux
    virtual void enableFlux()
    {
        computeFlux_ = true;
    }

    // Disable the computation of the flux
    virtual void disableFlux()
    {
        computeFlux_ = false;
        fluxPtr_.clear();
    }

    // Solve for the given unknown x, right-hand side source b and scalar
    // coefficient lambda the equation dx/dt = div(lambda*grad(x)) + b. Provide
    // lambda and b as pointers, so that they can be left empty which implies b
    // = 0 and lambda = 1. The ddt argument specifies whether or not to include
    // the d/dt term and the frac factor multiples the time step size.

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const meshField<Type,MeshType>* bPtr,
        const faceField<scalar,MeshType>* lambdaPtr,
        const bool ddt = false,
        const scalar dtFrac = 1
    ) = 0;

    // Single-argument support function

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const bool ddt = false,
        const scalar dtFrac = 1
    )
    {
        this->solve(x, nullptr, nullptr, ddt, dtFrac);
    }

    // Two-argument support function

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const meshField<Type,MeshType>* bPtr,
        const bool ddt = false,
        const scalar dtFrac = 1
    )
    {
        this->solve(x, bPtr, nullptr, ddt, dtFrac);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const meshField<Type,MeshType>& b,
        const bool ddt = false,
        const scalar dtFrac = 1
    )
    {
        this->solve(x, &b, nullptr, ddt, dtFrac);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const tmp<meshField<Type,MeshType>>& tb,
        const bool ddt = false,
        const scalar dtFrac = 1
    )
    {
        this->solve(x, &tb(), nullptr, ddt, dtFrac);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const faceField<scalar,MeshType>* lambdaPtr,
        const bool ddt = false,
        const scalar dtFrac = 1
    )
    {
        this->solve(x, nullptr, lambdaPtr, ddt, dtFrac);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const faceField<scalar,MeshType>& lambda,
        const bool ddt = false,
        const scalar dtFrac = 1
    )
    {
        this->solve(x, nullptr, &lambda, ddt, dtFrac);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const tmp<faceField<scalar,MeshType>>& tlambda,
        const bool ddt = false,
        const scalar dtFrac = 1
    )
    {
        this->solve(x, nullptr, &tlambda(), ddt, dtFrac);
    }

    // Three-argument support function

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const meshField<Type,MeshType>& b,
        const faceField<scalar,MeshType>* lambdaPtr,
        const bool ddt = false,
        const scalar dtFrac = 1
    )
    {
        this->solve(x, &b, lambdaPtr, ddt, dtFrac);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const tmp<meshField<Type,MeshType>>& tb,
        const faceField<scalar,MeshType>* lambdaPtr,
        const bool ddt = false,
        const scalar dtFrac = 1
    )
    {
        this->solve(x, &tb(), lambdaPtr, ddt, dtFrac);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const meshField<Type,MeshType>* bPtr,
        const faceField<scalar,MeshType>& lambda,
        const bool ddt = false,
        const scalar dtFrac = 1
    )
    {
        this->solve(x, bPtr, &lambda, ddt, dtFrac);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const meshField<Type,MeshType>& b,
        const faceField<scalar,MeshType>& lambda,
        const bool ddt = false,
        const scalar dtFrac = 1
    )
    {
        this->solve(x, &b, &lambda, ddt, dtFrac);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const tmp<meshField<Type,MeshType>>& tb,
        const faceField<scalar,MeshType>& lambda,
        const bool ddt = false,
        const scalar dtFrac = 1
    )
    {
        this->solve(x, &tb(), &lambda, ddt, dtFrac);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const meshField<Type,MeshType>* bPtr,
        const tmp<faceField<scalar,MeshType>>& tlambda,
        const bool ddt = false,
        const scalar dtFrac = 1
    )
    {
        this->solve(x, bPtr, &tlambda(), ddt, dtFrac);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const meshField<Type,MeshType>& b,
        const tmp<faceField<scalar,MeshType>>& tlambda,
        const bool ddt = false,
        const scalar dtFrac = 1
    )
    {
        this->solve(x, &b, &tlambda(), ddt, dtFrac);
    }

    virtual void solve
    (
        meshField<Type,MeshType>& x,
        const tmp<meshField<Type,MeshType>>& tb,
        const tmp<faceField<scalar,MeshType>>& tlambda,
        const bool ddt = false,
        const scalar dtFrac = 1
    )
    {
        this->solve(x, &tb(), &tlambda(), ddt, dtFrac);
    }
};

}

}

}

#ifdef NoRepository
    #include "PoissonSolver.C"
#endif

#endif
