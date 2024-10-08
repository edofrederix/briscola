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
        p.setOldTime();

        // Predictor

        USys = im::ddt(U);

        USys -= im::source(imSourceCoeff,U);
        USys -= exSource;

        USys -= im::laplacian(0.5*nu,U);
        USys -= ex::laplacian(0.5*nu,U);

        USys -= 0.5*H;
        H = ex::div(phi,U);
        USys += 1.5*H;

        // Solve predictor

        USolve->solve(USys + G);

        U += deltaT*G;
        U.correctBoundaryConditions();

        // Pressure equation

        phi = ex::faceFlux(U);

        Poisson->solve(p, ex::div(phi)/(-deltaT));

        G = ex::reconstruct(Poisson->flux());

        // Rhie-Chow correction

        U -= deltaT*G;
        U.correctBoundaryConditions();

        phi -= deltaT*Poisson->flux();

        io.write<colocated>();

        #include "colocatedContinuityErrors.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
}
