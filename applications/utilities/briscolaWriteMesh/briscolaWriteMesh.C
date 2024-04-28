#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"

#include "fvMesh.H"
#include "IO.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

int main(int argc, char *argv[])
{
    arguments::addBoolOption("ghosts", "Also write ghost cells");
    arguments::addBoolOption("partitioned", "Write data in partitioned format");

    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"
    #include "createBriscolaIO.H"

    runTime++;

    if (args.optionFound("ghosts"))
        io.enableGhosts();

    if (args.optionFound("partitioned"))
        io.enablePartitioned();

    forAll(fvMsh.msh(), l)
    {
        io.writeNow<colocated>(l);

        if (fvMsh.structured())
            io.writeNow<staggered>(l);
    }
}
