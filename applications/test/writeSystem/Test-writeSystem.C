#include "arguments.H"
#include "Time.H"
#include "fvMesh.H"
#include "linearSystem.H"
#include "imSchemes.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class MeshType>
void test(const fvMesh& fvMsh)
{
    meshField<scalar,MeshType> f
    (
        IOobject::groupName("f", MeshType::typeName),
        fvMsh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        true,
        true
    );

    f = Zero;
    f.correctBoundaryConditions();

    linearSystem<stencil,scalar,MeshType> sys1(im::laplacian(f));
    sys1.eliminateGhosts();

    OFstream file1(f.name() + "_laplacian");
    writeToFile(sys1,file1,0);

    linearSystem<diagStencil,scalar,MeshType> sys2(im::ddt(f));
    sys2.eliminateGhosts();

    OFstream file2(f.name() + "_ddt");
    writeToFile(sys2,file2,0);
}

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"

    test<colocated>(fvMsh);
    test<staggered>(fvMsh);
}
