#include "arguments.H"
#include "Time.H"
#include "uniformMesh.H"
#include "lineEdge.H"

using namespace Foam;
using namespace briscola;

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

    autoPtr<mesh> mshPtr(mesh::New(meshDict));

    const mesh& msh = mshPtr();

    // Find cell centers and see if they match

    forAll(msh, l)
    {
        const partLevel& lvl = msh[l];

        forAllBlock(lvl, i, j, k)
        {
            const vectorBlock points(lvl.points().cellPoints(i,j,k));

            const vector point(average(points));

            labelVector ijk = msh.findCell(point, l);

            if (ijk != labelVector(i,j,k))
            {
                FatalErrorInFunction
                    << "Test 1 failed" << endl << abort(FatalError);
            }

            if (!lvl.points().pointInCell(point, ijk))
                FatalErrorInFunction
                    << "Test 2 failed" << endl << abort(FatalError);
        }
    }
}
