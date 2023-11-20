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

    linearSystem<stencil,scalar,MeshType> sys(im::laplacian(f));
    sys.eliminateGhosts();

    OFstream file(f.name());

    writeToFile(sys,file,0);
}

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"

    test<colocated>(fvMsh);
    test<staggered>(fvMsh);
}
