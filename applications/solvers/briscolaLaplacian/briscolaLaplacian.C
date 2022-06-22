#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"

#include "fv.H"

#include <unistd.h>
#include <chrono>
using namespace std::chrono;

using namespace Foam;
using namespace briscola;
using namespace fv;

#define MESHTYPE colocated

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"

    // Solver dictionary

    IOdictionary solverDict
    (
        IOobject
        (
            "briscolaLaplacianDict",
            fvMsh.time().system(),
            fvMsh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    #include "createFields.H"
    #include "createBriscolaIO.H"

    while (runTime.run())
    {
        runTime++;

        Info << "Time = " << runTime.timeName() << endl;

        T.setOldTime();

        TSys = im::ddt(T);
        TSys -= im::laplacian(lambda,T);
        TSys -= source;

        TSolve->solve(TSys);

        io.write<MESHTYPE>();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
}
