#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"

#include "fv.H"
#include "vof.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

int main(int argc, char *argv[])
{
    arguments::addBoolOption("curvature", "Calculate interface curvature");

    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"

    Switch curvature = args.optionFound("curvature");

    #include "createBriscolaVof.H"

    #include "createFields.H"
    #include "createBriscolaIO.H"

    while (runTime.run())
    {
        #include "colocatedCourantNo.H"

        runTime++;

        Info << "Time = " << runTime.timeName() << endl;

        vf.solve(phi);

        if (curvature)
            curvatureSchemePtr->correct();

        io.write<colocated>();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
}
