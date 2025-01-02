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
void test(const fvMesh& fvMsh, const word solverType)
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

    OFstream os(f.name() + "_" + SType::typeName + "_" + solverType);
    sys.writeLevel(os);

    // Prepare solver and compute solution

    for (int nParts = 1; nParts <= Pstream::nProcs(); nParts++)
    {
        sys.x() = Zero;

        dictionary dict;

        dict.add("EigenSolver", solverType);
        dict.add("nAggregationParts", nParts);
        dict.add("maxIter", 100);
        dict.add("printStats", true);

        autoPtr<typename solver<SType,Type,MeshType>::directSolver> solverPtr
        (
            solver<SType,Type,MeshType>::directSolver::New
            (
                "Eigen",
                dict,
                fvMsh
            ).ptr()
        );

        solverPtr->prepare(sys);
        solverPtr->solve(sys);
    }

    // Write the solution

    OFstream oss
    (
        f.name()
      + "_"
      + SType::typeName
      + "_"
      + solverType
      + "_solution"
    );

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

    wordList types;

    types.append("default");
    types.append("partialPivLU");
    types.append("BiCGSTAB");

    #ifdef MKL
    types.append("Pardiso");
    #endif

    #ifdef SUITESPARSE
    types.append("UmfPack");
    #endif

    #ifdef SUPERLU
    types.append("SuperLU");
    #endif

    forAll(types, i)
    {
        const word type = types[i];

        test<symmStencil,scalar,colocated>(fvMsh, type);
        test<symmStencil,vector,colocated>(fvMsh, type);

        test<stencil,scalar,colocated>(fvMsh, type);
        test<stencil,vector,colocated>(fvMsh, type);

        test<stencil,scalar,staggered>(fvMsh, type);
        test<stencil,vector,staggered>(fvMsh, type);
    }
}
