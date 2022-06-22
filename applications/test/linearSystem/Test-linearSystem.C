#include "arguments.H"
#include "Time.H"

#include "linearSystem.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class SType, class Type, class MeshType>
void testConstructors(const fvMesh& fvMsh)
{
    meshField<Type,MeshType> f("f", fvMsh);
    meshField<Type,MeshType> g("g", fvMsh);

    // From a field

    linearSystem<SType,Type,MeshType> s1(f);

    // Copy constructors

    linearSystem<SType,Type,MeshType> s2a(s1);
    linearSystem<SType,Type,MeshType> s2b
    (
        tmp<linearSystem<SType,Type,MeshType>>
        (
            new linearSystem<SType,Type,MeshType>(f)
        )
    );

    // Copy constructor with other variable

    linearSystem<SType,Type,MeshType> s3a(s1, g);
    linearSystem<SType,Type,MeshType> s3b
    (
        tmp<linearSystem<SType,Type,MeshType>>
        (
            new linearSystem<SType,Type,MeshType>(f)
        ),
        g
    );
}

template<class SType, class Type, class MeshType>
void testResiduals(const fvMesh& fvMsh)
{
    meshField<Type,MeshType> f("f", fvMsh);
    meshField<Type,MeshType> res("res", fvMsh);

    linearSystem<SType,Type,MeshType> sys(f);

    sys.residual(res);
    res = sys.residual();

    forAll(f, l)
    {
        sys.residual(res[l]);
        res[l] = sys.residual(l);

        forAll(f[l], d)
        {
            sys.residual(res[l][d]);
            res[l][d] = sys.residual(l,d);
        }
    }
}

template<class SType, class Type, class MeshType>
void testMemberOperators(const fvMesh& fvMsh)
{
    meshField<Type,MeshType> f("f", fvMsh);
    meshField<Type,MeshType> g("g", fvMsh);

    linearSystem<SType,Type,MeshType> sys1(f);
    linearSystem<SType,Type,MeshType> sys2(f);

    sys1 = sys2;
    sys1 =
        tmp<linearSystem<SType,Type,MeshType>>
        (
            new linearSystem<SType,Type,MeshType>(f)
        );

    sys1 += sys2;
    sys1 +=
        tmp<linearSystem<SType,Type,MeshType>>
        (
            new linearSystem<SType,Type,MeshType>(f)
        );

    sys1 -= sys2;
    sys1 -=
        tmp<linearSystem<SType,Type,MeshType>>
        (
            new linearSystem<SType,Type,MeshType>(f)
        );

    sys1 = g;
    sys1 = (2*g);

    sys1 += g;
    sys1 += (2*g);

    sys1 -= g;
    sys1 -= (2*g);
}

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"

    IOdictionary meshDict
    (
        IOobject
        (
            runTime.system()/"briscolaMeshDict",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    fvMesh fvMsh(meshDict, runTime);

    // Colocated

    testConstructors<stencil,scalar,colocated>(fvMsh);
    testConstructors<stencil,vector,colocated>(fvMsh);
    testConstructors<symmStencil,scalar,colocated>(fvMsh);
    testConstructors<symmStencil,vector,colocated>(fvMsh);
    testConstructors<diagStencil,scalar,colocated>(fvMsh);
    testConstructors<diagStencil,vector,colocated>(fvMsh);

    testResiduals<stencil,scalar,colocated>(fvMsh);
    testResiduals<stencil,vector,colocated>(fvMsh);
    testResiduals<symmStencil,scalar,colocated>(fvMsh);
    testResiduals<symmStencil,vector,colocated>(fvMsh);
    testResiduals<diagStencil,scalar,colocated>(fvMsh);
    testResiduals<diagStencil,vector,colocated>(fvMsh);

    testMemberOperators<stencil,scalar,colocated>(fvMsh);
    testMemberOperators<stencil,vector,colocated>(fvMsh);
    testMemberOperators<symmStencil,scalar,colocated>(fvMsh);
    testMemberOperators<symmStencil,vector,colocated>(fvMsh);
    testMemberOperators<diagStencil,scalar,colocated>(fvMsh);
    testMemberOperators<diagStencil,vector,colocated>(fvMsh);

    // Staggered

    testConstructors<stencil,scalar,colocated>(fvMsh);
    testConstructors<stencil,vector,colocated>(fvMsh);
    testConstructors<symmStencil,scalar,colocated>(fvMsh);
    testConstructors<symmStencil,vector,colocated>(fvMsh);
    testConstructors<diagStencil,scalar,colocated>(fvMsh);
    testConstructors<diagStencil,vector,colocated>(fvMsh);

    testResiduals<stencil,scalar,colocated>(fvMsh);
    testResiduals<stencil,vector,colocated>(fvMsh);
    testResiduals<symmStencil,scalar,colocated>(fvMsh);
    testResiduals<symmStencil,vector,colocated>(fvMsh);
    testResiduals<diagStencil,scalar,colocated>(fvMsh);
    testResiduals<diagStencil,vector,colocated>(fvMsh);

    testMemberOperators<stencil,scalar,colocated>(fvMsh);
    testMemberOperators<stencil,vector,colocated>(fvMsh);
    testMemberOperators<symmStencil,scalar,colocated>(fvMsh);
    testMemberOperators<symmStencil,vector,colocated>(fvMsh);
    testMemberOperators<diagStencil,scalar,colocated>(fvMsh);
    testMemberOperators<diagStencil,vector,colocated>(fvMsh);
}
