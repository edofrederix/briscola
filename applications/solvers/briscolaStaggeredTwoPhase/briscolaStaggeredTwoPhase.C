#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"
#include "faceFluxScheme.H"
#include "interpolationScheme.H"

#include "fv.H"
#include "incompressibleTwoPhaseModel.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

int main(int argc, char *argv[])
{
    arguments::addBoolOption("FFT", "Use FFT for pressure equation");

    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"
    #include "createBriscolaTwoPhase.H"
    #include "createTimeControls.H"

    Switch split = args.optionFound("FFT");

    // This solver works for incompressible mixtures only

    incompressibleTwoPhaseModel& itpm =
        tpm.cast<incompressibleTwoPhaseModel>();

    #include "createRefs.H"
    #include "createFields.H"
    #include "createBriscolaIO.H"

    while (runTime.run())
    {
        #include "staggeredCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        const scalar deltaT = runTime.deltaTValue();
        const scalar deltaT0 = runTime.deltaT0Value();

        Info << "Time = " << runTime.timeName() << endl;

        U.setOldTime();

        // Update the volume fraction

        itpm.correct();

        // Predictor, Eq. (A.1) of Dodd & Ferrante (2014)

        USys = im::ddt(U);

        D = im::laplacian(mu,U)/rho;
        USys -= 0.5*D;
        USys -= 0.5*D.evaluate();

        USys -= 0.5*(deltaT/deltaT0)*H;

        phi = ex::faceFlux(U);

        H = ex::div(phi,U) - ((ex::grad(mu) & ex::grad(U)) / rho);

        USys += (1+0.5*(deltaT/deltaT0))*H;

        USys -= list(itpm.g());
        USys += ex::stagGrad(p) / rho - itpm.surfaceTension().stagForce();

        // Solve predictor

        USolve->solve(USys);

        // Pressure equation

        const staggeredScalarField rhoInv("rhoInv", 1.0/rho);
        const colocatedFaceScalarField coloRhoInvf("coloRhoInv",1.0/ex::interp(coloRho));

        if (split)
        {
            extrapolatedP = (1+deltaT/deltaT0)*p - (deltaT/deltaT0)*p.oldTime();
            extrapolatedP.correctBoundaryConditions();

            p.setOldTime();

            colocatedFaceScalarField splitCorrection = fa *
                (
                    (1 - coloMinRhof * coloRhoInvf) * ex::faceGrad(extrapolatedP)
                  + coloMinRhof * coloRhoInvf * ex::faceGrad(p)
                );

            Poisson->solve(p, coloMinRho*ex::coloDiv(U)/(-deltaT) - ex::div(splitCorrection));

            // Rhie-Chow correction

            U -= deltaT*
                (
                    ex::stagGrad(p)*minRhoInv
                  + (rhoInv - minRhoInv)*ex::stagGrad(extrapolatedP)
                  - ex::stagGrad(p.oldTime())*rhoInv
                );
            U.correctBoundaryConditions();
        }
        else
        {
            colocatedFaceScalarField splitCorrection = fa * coloRhoInvf * ex::faceGrad(p);
            p.setOldTime();

            Poisson->solve(p, ex::coloDiv(U)/(-deltaT) - ex::div(splitCorrection), coloRhoInvf);

            // Rhie-Chow correction

            U -= deltaT*(ex::stagGrad(p)-ex::stagGrad(p.oldTime()))*rhoInv;
            U.correctBoundaryConditions();
        }

        if (fvMsh.time().writeTime())
            Uc = ex::reconstruct(U);

        io.write<colocated>();
        io.write<staggered>();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
}
