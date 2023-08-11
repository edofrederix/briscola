#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"

#include "fv.H"
#include "immersedBoundary.H"

using namespace Foam;
using namespace briscola;
using namespace fv;
using namespace ibm;

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"
    #include "createTimeControls.H"

    // Solver dictionary

    IOdictionary solverDict
    (
        IOobject
        (
            "briscolaStaggeredDict",
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
        #include "staggeredCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        const scalar deltaT = runTime.deltaTValue();

        Info << "Time = " << runTime.timeName() << endl;

        U.setOldTime();

        // Predictor

        phi = ex::faceFlux(U);

        USys = im::ddt(U);
        USys += im::div(phi,U);
        USys -= im::laplacian(nu,U);
        USys += ex::stagGrad(p);
        USys -= source;

        // Immersed boundary

        if (solverDict.found("ImmersedBoundary"))
        {
            IB.penalization(USys);
            IB.IBM(USys);
        }

        // Solve momentum equation

        USolve->solve(USys);

        // Pressure equation

        U += deltaT*ex::stagGrad(p);
        U.correctBoundaryConditions();

        Poisson->solve(p, -ex::coloDiv(U)/deltaT);

        // Correction

        U -= deltaT*ex::stagGrad(p);
        U.correctBoundaryConditions();

        if (fvMsh.time().writeTime())
            Uc = ex::reconstruct(U);

        io.write<colocated>();
        io.write<staggered>();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
}
