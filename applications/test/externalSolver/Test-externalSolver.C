#include "arguments.H"
#include "Time.H"
#include "fvMesh.H"
#include "linearSystem.H"
#include "imSchemes.H"
#include "solver.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class SType, class Type, class MeshType>
void test(const fvMesh& fvMsh, const word solverType, const word subType)
{
    meshField<Type,MeshType> f
    (
        "f-"
      + word(MeshType::typeName)
      + "-"
      + word(pTraits<Type>::typeName),
        fvMsh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        true
    );

    f = 0.1*pTraits<Type>::one;

    linearSystem<SType,Type,MeshType> sys(im::laplacian<SType>(f));
    sys -= im::ddt(f);

    sys.eliminateGhosts();
    sys.singular();
    sys.diagonal();

    // Write the system to a file

    const word fileName
    (
        f.name() + "_" +
        SType::typeName + "_" +
        solverType + "_" +
        subType
    );

    Info<< "Writing to " << fileName << endl;

    OFstream os(fileName);
    sys.writeLevel(os);

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

        autoPtr<typename solver<SType,Type,MeshType>::externalSolver> solverPtr
        (
            solver<SType,Type,MeshType>::externalSolver::New
            (
                solverType,
                dict,
                fvMsh
            ).ptr()
        );

        solverPtr->prepare(sys);
        solverPtr->solve(sys);
    }

    // Write last solution

    OFstream oss(word(fileName + "_solution"));

    for (int d = 0; d < MeshType::numberOfDirections; d++)
    {
        List<List<Type>> data(Pstream::nProcs());
        data[Pstream::myProcNo()].setSize(cmptProduct(fvMsh.N<MeshType>(d)));

        int l = 0;
        forAllCells(sys.x().direction(d), i, j, k)
            data[Pstream::myProcNo()][l++] =
                sys.x().direction(d)(i,j,k);

        Pstream::gatherList(data);

        if (Pstream::master())
            forAll(data, proc)
                forAll(data[proc], j)
                    for (int i = 0; i < pTraits<Type>::nComponents; i++)
                        oss << scalar_cast(&data[proc][j])[i]
                            << (i == pTraits<Type>::nComponents-1 ? nl : ' ');
    }
}

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
    subTypes[findIndex(solverTypes,"Eigen")].append("PartialPivLU");
    subTypes[findIndex(solverTypes,"Eigen")].append("BiCGSTAB");

    #ifdef SUPERLU
    subTypes[findIndex(solverTypes,"PETSc")].append("SuperLU");
    subTypes[findIndex(solverTypes,"Eigen")].append("SuperLU");
    #endif

    #ifdef SUPERLU_DIST
    subTypes[findIndex(solverTypes,"PETSc")].append("SuperLUDist");
    #endif

    forAll(solverTypes, i)
    {
        const word solverType = solverTypes[i];

        forAll(subTypes[i], j)
        {
            const word subType = subTypes[i][j];

            test<stencil,scalar,colocated>(fvMsh, solverType, subType);
            test<stencil,vector,colocated>(fvMsh, solverType, subType);

            test<stencil,scalar,staggered>(fvMsh, solverType, subType);
            test<stencil,vector,staggered>(fvMsh, solverType, subType);
        }
    }
}
