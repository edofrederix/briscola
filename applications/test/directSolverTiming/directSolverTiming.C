#include "arguments.H"
#include "Time.H"
#include "fvMesh.H"
#include "linearSystem.H"
#include "imSchemes.H"
#include "solver.H"
#include "ticToc.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

typedef typename solver<stencil,scalar,colocated>::directSolver
    directSolverType;

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"

    wordList types;

    types.append("SparseLU");
    types.append("BiCGSTAB");

    #ifdef SUPERLU
    types.append("SuperLU");
    #endif

    colocatedScalarField f
    (
        "f",
        fvMsh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        true
    );

    f = 0.1;

    colocatedScalarSystem sys(im::laplacian(f));
    sys -= im::ddt(f);

    sys.x().restrict();
    sys.b().restrict();

    sys.eliminateGhosts();
    sys.singular();
    sys.diagonal();

    forAll(types, i)
    {
        const word type = types[i];

        // Prepare solver and compute solution

        for (int nParts = 1; nParts <= Pstream::nProcs(); nParts++)
        {
            sys.x() = Zero;

            dictionary dict;

            dict.add("EigenSolver", type);
            dict.add("nAggregationParts", nParts);
            dict.add("maxIter", 100);
            dict.add("printStats", true);

            forAll(sys.x(), l)
            {
                Info<< "Solver = " << type
                    << ", nParts = " << nParts
                    << ", level = " << l << endl;

                initTicToc(2)

                autoPtr<directSolverType>
                    solverPtr
                    (
                        directSolverType::New
                        (
                            "Eigen",
                            dict,
                            fvMsh,
                            l
                        ).ptr()
                    );

                tic(0)
                solverPtr->prepare(sys);
                toc(0)

                tic(1)
                solverPtr->solve(sys);
                toc(1)

                printTicToc
            }
        }
    }
}
