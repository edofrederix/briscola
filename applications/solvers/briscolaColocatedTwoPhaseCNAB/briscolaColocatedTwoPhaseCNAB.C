#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"

#include "fv.H"
#include "TwoPhaseModel.H"

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
    #include "createBriscolaColocatedTwoPhase.H"
    #include "createTimeControls.H"

    #include "createRefs.H"
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

        // Update the two-phase model and specific volumes

        twoPhase.correct();

        v = 1.0/rho;
        vf = ex::interp(v);

        // Predictor

        USys = im::ddt(U);

        USys -= im::source(imSourceCoeff,U)*v;
        USys -= exSource*v;

        USys -= 0.5*H;
        H = ex::div(phi,U);
        USys += 1.5*H;

        if (modified)
        {
            USys -= im::laplacian(9.0/16.0*mu,U)*v;

            USys -= 1.0/16.0*G;
            G = ex::laplacian(mu,U)*v;
            USys -= 6.0/16.0*G;
        }
        else
        {
            USys -= im::laplacian(0.5*mu,U)*v;
            USys -= ex::laplacian(0.5*mu,U)*v;
        }

        USys -= twoPhase.buoyancy()*v;

        // Solve predictor

        USolve->solve(USys);

        // Pressure equation

        colocatedFaceScalarField phiStar
        (
            ex::faceFlux(U)
          + deltaT*twoPhase.flux()*vf
        );

        Poisson->solve(p, ex::div(phiStar)/(-deltaT), vf);

        // Rhie-Chow correction

        U -=
            deltaT
          * ex::reconstruct
            (
                Poisson->flux()/vf
              - twoPhase.flux()
            )*v;

        U.correctBoundaryConditions();

        phi = phiStar - deltaT*Poisson->flux();

        io.write<colocated>();

        #include "colocatedContinuityErrors.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
}
