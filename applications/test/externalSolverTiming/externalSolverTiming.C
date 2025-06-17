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

typedef typename solver<stencil,scalar,colocated>::externalSolver
    externalSolverType;

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"

    wordList solverTypes;

    solverTypes.append("PETSc");
    solverTypes.append("Eigen");

    List<wordList> subTypes(solverTypes.size());

    subTypes[findIndex(solverTypes,"PETSc")].append("PCLU");
    subTypes[findIndex(solverTypes,"PETSc")].append("KSPBCGS");
    subTypes[findIndex(solverTypes,"PETSc")].append("KSPIBCGS");
    subTypes[findIndex(solverTypes,"PETSc")].append("KSPGMRES");
    subTypes[findIndex(solverTypes,"PETSc")].append("KSPFGMRES");

    subTypes[findIndex(solverTypes,"Eigen")].append("SparseLU");
    subTypes[findIndex(solverTypes,"Eigen")].append("BiCGSTAB");

    #ifdef SUPERLU
    subTypes[findIndex(solverTypes,"PETSc")].append("SuperLU");
    subTypes[findIndex(solverTypes,"Eigen")].append("SuperLU");
    #endif

    #ifdef SUPERLU_DIST
    subTypes[findIndex(solverTypes,"PETSc")].append("SuperLUDist");
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

    forAll(solverTypes, i)
    {
        const word solverType = solverTypes[i];

        forAll(subTypes[i], j)
        {
            const word subType = subTypes[i][j];

            // Prepare solver and compute solution

            for (int nParts = 1; nParts <= Pstream::nProcs(); nParts++)
            {
                sys.x() = Zero;

                dictionary dict;

                dict.add(word(solverType + "Solver"), subType);

                dict.add("maxIter", 100);
                dict.add("relTol", 1e-12);
                dict.add("tolerance", 1e-12);
                dict.add("printStats", false);
                dict.add("nAggregationParts", nParts);

                forAll(sys.x(), l)
                {
                    Info<< "Solver = " << solverType
                        << ", sub-solver = " << subType
                        << ", nParts = " << nParts
                        << ", level = " << l << " ";

                    initTicToc(3)

                    autoPtr<externalSolverType>
                        solverPtr
                        (
                            externalSolverType::New
                            (
                                solverType,
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

                    tic(2)
                    solverPtr->solve(sys);
                    toc(2)

                    printTicToc
                }
            }
        }
    }
}
