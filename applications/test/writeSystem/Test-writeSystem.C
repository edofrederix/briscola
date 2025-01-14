#include "arguments.H"
#include "Time.H"
#include "fvMesh.H"
#include "linearSystem.H"
#include "imSchemes.H"
#include "restrictionScheme.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class SType, class Type, class MeshType>
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
        true
    );

    f = 0.1*pTraits<Type>::one;

    linearSystem<SType,Type,MeshType> sys(im::laplacian<SType>(f));
    sys -= im::ddt(f);

    sys.eliminateGhosts();
    sys.singular();
    sys.diagonal();

    sys.x().restrict();
    sys.b().restrict();

    // Write all levels

    forAll(f, l)
    {
        OFstream os(f.name() + "_" + Foam::name(l));
        sys.writeLevel(os, l);
    }
}

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"

    test<stencil,scalar,colocated>(fvMsh);
    test<stencil,vector,colocated>(fvMsh);

    test<stencil,scalar,staggered>(fvMsh);
    test<stencil,vector,staggered>(fvMsh);
}
