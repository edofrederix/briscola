#include "arguments.H"
#include "Time.H"
#include "fvMesh.H"
#include "linearSystem.H"
#include "imSchemes.H"
#include "restrictionScheme.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class SType, class MeshType>
void test(const fvMesh& fvMsh)
{
    meshField<scalar,MeshType> f
    (
        word("f-" + word(MeshType::typeName)),
        fvMsh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        true
    );

    f = 0.1;

    linearSystem<SType,scalar,MeshType> sys(im::laplacian<SType>(f));
    sys -= im::ddt(f);
    sys.eliminateGhosts();

    // Write

    writeToFile(sys, f.name(), 0);
}

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"

    test<symmStencil,colocated>(fvMsh);
    test<stencil,staggered>(fvMsh);
}
