#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"
#include "immersedBoundaryConditionStaggeredMassSource.C"

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

        Info << "Time = " << runTime.name() << endl;

        U.setOldTime();
        p.setOldTime();

        phi = ex::faceFlux(U);

        // Predictor

        USys = im::ddt(U);

        USys -= im::source(imSourceCoeff,U);
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

        Poisson->solve(p, ibmCorr(ex::coloDiv(U),U)/(-deltaT));

        // Correction

        U -= deltaT*ex::stagGrad(p);
        U.correctBoundaryConditions();

        // Reconstruct the colocated velocity

        Uc = ex::reconstruct(U);
        Uc.correctBoundaryConditions();

        io.write<colocated>();
        io.write<staggered>();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
}
