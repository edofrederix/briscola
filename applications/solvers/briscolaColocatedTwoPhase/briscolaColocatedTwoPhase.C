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
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"
    #include "createBriscolaTwoPhase.H"
    #include "createTimeControls.H"

    // This solver works for incompressible mixtures only

    incompressibleTwoPhaseModel& itpm =
        tpm.cast<incompressibleTwoPhaseModel>();

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

        // Update the volume fraction

        itpm.correct();

        // Predictor, Eq. (A.1) of Dodd & Ferrante (2014)

        USys = im::ddt(U);

        D = im::laplacian(mu,U)/rho;
        USys -= 0.5*D;
        USys -= 0.5*D.evaluate();

        USys -= 0.5*H;
        H = (ex::div(phi,U) - itpm.surfaceTension());
        USys += 1.5*H;

        USys -= itpm.g();

        // Solve predictor

        USolve->solve(USys);

        // Pressure equation

        phi = ex::faceFlux(U);

        const colocatedScalarField rhoInv("rhoInv", 1.0/rho);
        const colocatedLowerFaceScalarField rhoInvf
        (
            "rhoInv",
            1.0/ex::interp(rho)
        );

        Poisson->solve(p, ex::div(phi)/(-deltaT), rhoInvf);

        // Rhie-Chow correction

        U -= deltaT*ex::grad(p)*rhoInv;
        U.correctBoundaryConditions();

        phi -= deltaT*ex::faceGrad(p)*fa*rhoInvf;

        io.write<colocated>();

        #include "colocatedContinuityErrors.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
}
