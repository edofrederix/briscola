#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"

#include "fv.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

// Tool to test restriction by setting cell-centered coordinates at the finest
// level, restricting those and writing all levels. The outcome of the selected
// restriction operator (defaulting to linear) can then be viewed.

template<class Type, class MeshType>
void restrict(meshField<Type,MeshType>& field, const word scheme)
{
    autoPtr<restrictionScheme<vector,MeshType>> schemePtr
    (
        restrictionScheme<vector,MeshType>::New
        (
            field.fvMsh(),
            scheme
        )
    );

    for (int l = 0; l < field.fvMsh().size()-1; l++)
    {
        schemePtr->restrict(field[l+1], field[l]);
        field[l+1].correctBoundaryConditions();
    }
}

int main(int argc, char *argv[])
{
    arguments::addOption
    (
        "scheme",
        "restriction scheme",
        "specify the restriction scheme (defaults to linear)"
    );

    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"
    #include "createBriscolaIO.H"

    runTime++;

    word scheme =
        args.optionFound("scheme")
      ? word(args["scheme"])
      : word("linear");

    colocatedVectorField C
    (
        "C",
        fvMsh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true,
        true
    );

    C = Zero;

    C[0] = fvMsh.metrics<colocated>().cellCenters()[0];

    restrict(C, scheme);

    autoPtr<meshField<vector,staggered>> SPtr;

    if (fvMsh.structured())
    {
        SPtr.reset
        (
            new staggeredVectorField
            (
                "S",
                fvMsh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                true,
                true
            )
        );

        meshField<vector,staggered>& S = SPtr();

        S = Zero;

        S[0] = fvMsh.metrics<staggered>().cellCenters()[0];

        restrict(S, scheme);
    }

    forAll(fvMsh.msh(), l)
    {
        io.writeNow<colocated>(l);
        io.writeNow<staggered>(l);
    }
}
