#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"

#include "fv.H"
#include "incompressibleTwoPhaseModel.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

int main(int argc, char *argv[])
{
    arguments::addBoolOption("split", "Split the pressure equation");

    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"
    #include "createBriscolaTwoPhase.H"
    #include "createTimeControls.H"

    // This solver works for incompressible mixtures only

    incompressibleTwoPhaseModel& icoTwoPhase =
        twoPhase.cast<incompressibleTwoPhaseModel>();

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
        const scalar deltaT0 = runTime.deltaT0Value();

        Info << "Time = " << runTime.timeName() << endl;

        U.setOldTime();
        p.setOldTime();

        // Update the two-phase model and specific volumes

        icoTwoPhase.correct();

        v = 1.0/rho;
        vf = ex::interp(v);

        // Predictor

        USys = im::ddt(U);

        USys -= im::source(imSourceCoeff,U);
        USys -= exSource;

        USys -= im::laplacian(mu,U,0.5)*v;
        USys -= ex::laplacian(mu,U,0.5)*v;

        USys -= 0.5*(deltaT/deltaT0)*H;

        H = ex::div(phi,U)
          - ex::div(mu*ex::faceFlux(dev2(T(ex::grad(U)))))*v
          - ex::reconstruct(icoTwoPhase.surfaceTension())*v;

        USys += (1.0 + 0.5*(deltaT/deltaT0))*H;
        USys -= icoTwoPhase.g();

        // Solve predictor

        USolve->solve(USys + ex::grad(p)*v);

        U += deltaT*ex::grad(p)*v;
        U.correctBoundaryConditions();

        // Pressure equation

        phi = ex::faceFlux(U);

        Poisson->solve(p, ex::div(phi)/(-deltaT), vf);

        // Rhie-Chow correction

        U -= deltaT*ex::reconstruct(Poisson->flux()/vf)*v;
        U.correctBoundaryConditions();

        phi -= deltaT*Poisson->flux();

        io.write<colocated>();

        #include "colocatedContinuityErrors.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
}
