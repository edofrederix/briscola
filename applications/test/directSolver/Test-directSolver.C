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
        true,
        true
    );

    f = 0.1*pTraits<Type>::one;

    autoPtr<typename solver<SType,Type,MeshType>::directSolver> solverPtr
    (
        solver<SType,Type,MeshType>::directSolver::New
        (
            solverType,
            dictionary::null,
            fvMsh
        ).ptr()
    );

    linearSystem<SType,Type,MeshType> sys(im::laplacian<SType>(f));
    sys -= im::ddt(f);

    sys.eliminateGhosts();
    sys.singular();
    sys.diagonal();

    // Write the system to a file

    writeToFile(sys, f.name() + "_" + SType::typeName + "_" + solverType);

    // Prepare solver and compute solution

    solverPtr->prepare(sys);
    solverPtr->solve(sys);

    // Write the solution

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
        {
            const fileName name =
                f.name()
              + "_"
              + SType::typeName
              + "_"
              + solverType
              + (
                    MeshType::numberOfDirections > 1
                  ? "_" + Foam::name(d)
                  : ""
                )
              + "_solution";

            OFstream file(name);

            forAll(data, proc)
                forAll(data[proc], l)
                    file<< data[proc][l] << nl;
        }
    }
}

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"

    wordList types(2);

    types[0] = "APLU";
    types[1] = "Eigen";

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
