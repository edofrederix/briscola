#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"

#include "fv.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

// Crank-Nicolson/Adams-Bashforth (CNAB) IMEX scheme, see Eq. (2) of Ascher, U.
// M., Ruuth, S. J., & Wetton, B. T. (1995). Implicit-explicit methods for
// time-dependent partial differential equations. SIAM Journal on Numerical
// Analysis, 32(3), 797-823.

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

        Info << "Time = " << runTime.name() << endl;

        U.setOldTime();
        p.setOldTime();

        // Predictor

        USys = im::ddt(U);

        USys -= exSource;

        USys -= 0.5*H;
        H = ex::div(phi,U);
        USys += 1.5*H;

        if (modified)
        {
            USys -= im::laplacian(9.0/16.0*nu,U);

            USys -= 1.0/16.0*G;
            G = ex::laplacian(nu,U);
            USys -= 6.0/16.0*G;
        }
        else
        {
            USys -= im::laplacian(0.5*nu,U);
            USys -= ex::laplacian(0.5*nu,U);
        }

        // Solve predictor

        USolve->solve(USys);

        // Pressure equation

        phi = ex::faceFlux(U);

        Poisson->solve(p, ex::div(phi)/(-deltaT));

        // Rhie-Chow correction

        U -= deltaT*ex::grad(p);
        U.correctBoundaryConditions();

        phi -= deltaT*Poisson->flux();

        io.write<colocated>();

        #include "colocatedContinuityErrors.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
}
