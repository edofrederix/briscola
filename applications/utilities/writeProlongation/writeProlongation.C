#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"

#include "fv.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

// Tool to test prolongation by setting cell-centered coordinates at the
// coarsest level, prolonging those and writing all levels. The outcome of the
// selected prolongation operator (defaulting to linear) can then be viewed.

template<class Type, class MeshType>
void prolong(meshField<Type,MeshType>& field, const word scheme)
{
    autoPtr<prolongationScheme<vector,MeshType>> schemePtr
    (
        prolongationScheme<vector,MeshType>::New
        (
            field.fvMsh(),
            scheme
        )
    );

    for (int l = field.fvMsh().size()-1; l >= 1; l--)
        for (int d = 0; d < MeshType::numberOfDirections; d++)
            schemePtr->prolong(field[l-1][d], field[l][d], eqOp<Type>());
}

int main(int argc, char *argv[])
{
    arguments::addOption
    (
        "scheme",
        "prolongation scheme",
        "specify the prolongation scheme (defaults to linear)"
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

    const label nLevels = fvMsh.size();

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

    C[nLevels-1] =
        fvMsh.metrics<colocated>().cellCenters()[nLevels-1];

    prolong(C, scheme);

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

        S[nLevels-1] =
            fvMsh.metrics<staggered>().cellCenters()[nLevels-1];

        prolong(S, scheme);
    }

    forAll(fvMsh.msh(), l)
    {
        io.writeNow<colocated>(l);
        io.writeNow<staggered>(l);
    }
}
