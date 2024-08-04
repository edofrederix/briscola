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

        icoTwoPhase.correct();

        v = 1.0/rho;
        vcf = ex::coloFaceInterp(v);

        // Predictor, Eq. (A.1) of Dodd & Ferrante (2014)

        USys = im::ddt(U);

        USys -= im::source(imSourceCoeff,U);
        USys -= v*exSource;

        USys -= im::laplacian(mu,U,0.5)*v;
        USys -= ex::laplacian(mu,U,0.5)*v;

        USys -= 0.5*(deltaT/deltaT0)*H;

        phi = ex::faceFlux(U);

        H = ex::div(phi,U)
          - stagDotProduct(ex::grad(mu),ex::grad(U))*v
          - ex::stagReconstruct(icoTwoPhase.surfaceTension())*v;

        USys += (1.0 + 0.5*(deltaT/deltaT0))*H;

        if (reduced)
        {
            USys -= ex::stagReconstruct(icoTwoPhase.buoyancy())*v;
        }
        else
        {
            USys -= list(icoTwoPhase.g());
        }

        if (fvMsh.immersedBoundaryPresent())
        {
            U.correctImmersedBoundaryConditions(USys);
        }

        for (int corr = 0; corr < nCorr; corr++)
        {
            // Solve predictor with latest pressure

            USolve->solve(USys + ex::stagGrad(p)*v);

            USysAux = USys + ex::stagGrad(p)*v;

            residual = USysAux.residual();

            // Pressure equation

            if (split)
            {
                q = (1.0 + deltaT/deltaT0)*p - (deltaT/deltaT0)*p.oldTime();
                p.setOldTime();

                colocatedLowerFaceScalarField corr
                (
                    fa
                  * (
                        ex::faceGrad(q)*(1.0 - minRho*vcf)
                      + ex::faceGrad(p)*minRho*vcf
                    )
                );

                Poisson->solve
                (
                    p,
                    minRho*ex::coloDiv(U)/(-deltaT) - ex::div(corr)
                );

                // Correct velocity

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

                Poisson->solve
                (
                    p,
                    ex::coloDiv(U)/(-deltaT) - ex::div(corr),
                    vcf
                );

                // Correct velocity

                U -= deltaT*(ex::stagGrad(p) - ex::stagGrad(p.oldTime()))*v;
                U.correctBoundaryConditions();
            }
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
