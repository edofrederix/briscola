#include "arguments.H"
#include "Time.H"
#include "fvMesh.H"
#include "linearSystem.H"
#include "imSchemes.H"
#include "directSolver.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class Type, class MeshType>
void test(const fvMesh& fvMsh)
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

    dictionary dict;
    dict.add("type", "APLU");

    autoPtr<directSolver<stencil,Type,MeshType>> solverPtr
    (
        directSolver<stencil,Type,MeshType>::New(dict,fvMsh).ptr()
    );

    linearSystem<stencil,Type,MeshType> sys(im::laplacian(f));
    sys -= im::ddt(f);
    sys.eliminateGhosts();

    // Write the system to a file

    writeToFile(sys,f.name());

    // Compute the solution

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

    test<scalar,colocated>(fvMsh);
    test<scalar,staggered>(fvMsh);

    test<vector,colocated>(fvMsh);
    test<vector,staggered>(fvMsh);
}
