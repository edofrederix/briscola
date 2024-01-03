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
    arguments::addBoolOption("split", "Split the pressure equation");

    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"
    #include "createBriscolaTwoPhase.H"
    #include "createTimeControls.H"

    Switch split = args.optionFound("split");

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

        // Update the two-phase model and specific volumes

        itpm.correct();

        v = itpm.v<staggered>();
        vc = itpm.v<colocated>();
        vcf = ex::interp(vc);

        // Predictor, Eq. (A.1) of Dodd & Ferrante (2014)

        USys = im::ddt(U);

        D = im::laplacian<stencil>(mu,U)/rho;
        USys -= 0.5*D;
        USys -= 0.5*D.evaluate();

        USys -= 0.5*(deltaT/deltaT0)*H;

        phi = ex::faceFlux(U);

        gradMu = ex::grad(mu);
        gradU = ex::grad(U);

        H = ex::div(phi,U) - ex::transposeMultiplicate(gradMu,gradU)/rho;

        USys += (1.0 + 0.5*(deltaT/deltaT0))*H;

        USys -= list(itpm.g());
        USys += ex::stagGrad(p)/rho - itpm.surfaceTension().stagForce();

        // Solve predictor

        USolve->solve(USys);

        // Pressure equation

        if (split)
        {
            q = (1.0 + deltaT/deltaT0)*p - (deltaT/deltaT0)*p.oldTime();
            p.setOldTime();

            colocatedLowerFaceScalarField corr
            (
                fa
              * (
                    ex::faceGrad(q) * (1.0 - minRho*vcf)
                  + ex::faceGrad(p) * minRho*vcf
                )
            );

            Poisson->solve(p, minRho*ex::coloDiv(U)/(-deltaT) - ex::div(corr));

            // Rhie-Chow correction

            U -=
                deltaT
              * (
                    ex::stagGrad(p)*maxv
                  - ex::stagGrad(p.oldTime())*v
                  + ex::stagGrad(q)*(v - maxv)
                );

            U.correctBoundaryConditions();
        }
        else
        {
            colocatedLowerFaceScalarField corr(fa*vcf*ex::faceGrad(p));
            p.setOldTime();

            Poisson->solve(p, ex::coloDiv(U)/(-deltaT) - ex::div(corr), vcf);

            // Rhie-Chow correction

            U -= deltaT*(ex::stagGrad(p) - ex::stagGrad(p.oldTime()))*v;
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
