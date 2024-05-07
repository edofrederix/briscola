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

        // Update the two-phase model and specific volumes

        icoTwoPhase.correct();

        v = 1.0/rho;
        vf = ex::interp(v);

        // Predictor, Eq. (A.1) of Dodd & Ferrante (2014)

        USys = im::ddt(U);

        USys -= im::source(imSourceCoeff,U);
        USys -= exSource;

        USys -= im::laplacian(mu,U,0.5)*v;
        USys -= ex::laplacian(mu,U,0.5)*v;

        USys -= 0.5*(deltaT/deltaT0)*H;

        H = ex::div(phi,U)
          - (ex::grad(mu) & ex::grad(U))*v
          - ex::reconstruct(icoTwoPhase.surfaceTension())*v;

        USys += (1.0 + 0.5*(deltaT/deltaT0))*H;

        if (reduced)
        {
            USys -= ex::reconstruct(icoTwoPhase.buoyancy())*v;
        }
        else
        {
            USys -= icoTwoPhase.g();
        }

        for (int corr = 0; corr < nCorr; corr++)
        {
            // Solve predictor

            USolve->solve(USys + ex::grad(p)*v);

            // Pressure equation

            phi = ex::faceFlux(U);

            if (split)
            {
                q = (1.0 + deltaT/deltaT0)*p - deltaT/deltaT0*p.oldTime();
                p.setOldTime();

                colocatedLowerFaceScalarField corr
                (
                    fa
                  * (
                        ex::faceGrad(q)*(1.0 - minRho*vf)
                      + ex::faceGrad(p)*minRho*vf
                    )
                );

                Poisson->solve
                (
                    p,
                    minRho*ex::div(phi)/(-deltaT) - ex::div(corr)
                );

                // Rhie-Chow correction

                U -=
                    deltaT
                  * (
                        ex::grad(p)*maxv
                      - ex::grad(p.oldTime())*v
                      + ex::grad(q)*(v - maxv)
                    );

                U.correctBoundaryConditions();

                phi -=
                    deltaT*fa
                  * (
                        ex::faceGrad(p)*maxv
                      - ex::faceGrad(p.oldTime())*vf
                      + ex::faceGrad(q)*(vf - maxv)
                    );
            }
            else
            {
                colocatedLowerFaceScalarField corr(fa*vf*ex::faceGrad(p));
                p.setOldTime();

                Poisson->solve(p, ex::div(phi)/(-deltaT) - ex::div(corr), vf);

                // Rhie-Chow correction

                U -= deltaT*(ex::grad(p) - ex::grad(p.oldTime()))*v;
                U.correctBoundaryConditions();

                phi -=
                    deltaT*(ex::faceGrad(p) - ex::faceGrad(p.oldTime()))*fa*vf;
            }
        }

        io.write<colocated>();

        #include "colocatedContinuityErrors.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
}
