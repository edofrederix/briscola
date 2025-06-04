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
    #include "createBriscolaStaggeredTwoPhase.H"
    #include "createTimeControls.H"

    #include "createRefs.H"
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
        p.setOldTime();

        phi = ex::faceFlux(U);

        // Update the two-phase model and specific volumes

        twoPhase.correct();

        v = 1.0/rho;
        vcf = ex::coloFaceInterp(v);

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

        U += deltaT*ex::stagReconstruct(twoPhase.flux())*v;
        U.correctBoundaryConditions();

        // Pressure equation

        Poisson->solve(p, ex::coloDiv(U)/(-deltaT), vcf);

        // Correct velocity

        U -= deltaT*ex::stagReconstruct(Poisson->flux()/vcf)*v;
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
