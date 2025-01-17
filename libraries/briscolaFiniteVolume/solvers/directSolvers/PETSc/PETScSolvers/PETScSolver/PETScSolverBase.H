#ifndef PETScSolverBase_H
#define PETScSolverBase_H

#include "dictionary.H"
#include "autoPtr.H"
#include <petscksp.h>

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class PETScLinearSystem;

class PETScSolverBase
{
protected:

    // Copy of the dictionary
    dictionary dict_;

public:

    // Runtime type information
    TypeName("PETScSolverBase");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        PETScSolverBase,
        dictionary,
        (const dictionary& dict),
        (dict)
    );

    // Construct from dictionary
    PETScSolverBase(const dictionary& dict);

    // Copy constructor
    PETScSolverBase(const PETScSolverBase&);

    // Selector from solver type and coefficient dictionary
    static autoPtr<PETScSolverBase> New(const word, const dictionary& dict);

    // Destructor
    virtual ~PETScSolverBase();

    // Prepare from a sparse matrix
    virtual void prepare(const Mat& mat, const MPI_Comm& comm) = 0;

    // Solve given a right-hand side
    virtual void solve(Vec& x, const Vec& rhs, const MPI_Comm& comm) = 0;

    // Parallel solver bool
    virtual bool parallel() const = 0;
};

}

}

}

#endif
