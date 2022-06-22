#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"

#include "fv.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

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
            "briscolaColocatedDict",
            fvMsh.time().system(),
            fvMsh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    #include "createFields.H"
    #include "createRefs.H"
    #include "createBriscolaIO.H"
    #include "initContinuityErrors.H"

    while (runTime.run())
    {
        #include "colocatedCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        const scalar deltaT = runTime.deltaTValue();

        Info << "Time = " << runTime.timeName() << endl;

        U.setOldTime();

        // Predictor

        USys = im::ddt(U);
        USys += im::div(phi,U);
        USys -= im::laplacian(nu,U);
        USys += ex::grad(p);
        USys -= source;

        USolve->solve(USys);

        // Pressure equation

        U += deltaT*ex::grad(p);
        U.correctBoundaryConditions();

        phi = ex::faceFlux(U);

        pSys = im::laplacian(p);
        pSys -= ex::div(phi)/deltaT;

        pSolve->solve(pSys);

        // Rhie-Chow correction

        U -= deltaT*ex::grad(p);
        U.correctBoundaryConditions();

        phi -= deltaT*ex::faceGrad(p)*fa;

        io.write<colocated>();

        #include "colocatedContinuityErrors.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
}
