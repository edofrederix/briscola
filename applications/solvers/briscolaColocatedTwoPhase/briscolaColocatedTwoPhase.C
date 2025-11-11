#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"

#include "fv.H"
#include "TwoPhaseModel.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"
    #include "createRungeKuttaScheme.H"
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

        Info << "Time = " << runTime.name() << endl;

        U.setOldTime();
        p.setOldTime();

        // Update the two-phase model and specific volumes

        twoPhase.correct();

        v = 1.0/rho;
        vf = ex::interp(v);
        vf.max(1e-12);

        // Explicit source

        tmp<colocatedVectorField> tSource =
            (exSource + twoPhase.buoyancy())*v;

        while (rk.loop())
        {
            const label stage = rk.stage();

            if (rk.solve())
            {
                const scalar A = rk.A();
                const scalar B = rk.B();
                const scalar C = rk.C();

                // Predictor

                USys = im::ddt(U);
                USys -= C*tSource();
                USys -= rk.stageSum(stageSourcesA, stageSourcesB);

                if (rk.imStageA())
                {
                    tUSysA = -im::div(phi,U);
                    USys -= A*tUSysA.ref();
                }

                if (rk.imStageB())
                {
                    tUSysB =
                        (
                            im::laplacian(mu,U)
                          + ex::div(mu*ex::faceDotGrad(U))
                        )*v;

                    USys -= B*tUSysB.ref();
                }

                // Solve predictor

                USolve->solve(USys);

                // Pressure equation

                colocatedScalarFaceField phiStar
                (
                    ex::faceFlux(U)
                  + C*deltaT*twoPhase.flux()*vf
                );

                Poisson->solve(p, ex::div(phiStar)/(-C*deltaT), vf);

                // Rhie-Chow correction

                U -=
                    C*deltaT
                  * ex::reconstruct
                    (
                        Poisson->flux()/vf
                      - twoPhase.flux()
                    )*v;

                U.correctBoundaryConditions();

                if (rk.lastStage())
                    phi = phiStar - C*deltaT*Poisson->flux();
            }

            // Store Runge-Kutta sources

            if (rk.storeStageA())
                stageSourcesA[stage-1] =
                    rk.solve() && rk.imStageA()
                  ? tUSysA->evaluate()
                  : -ex::div(phi,U);

            if (rk.storeStageB())
                stageSourcesB[stage-1] =
                    rk.solve() && rk.imStageB()
                  ? tUSysB->evaluate()
                  : v*(ex::laplacian(mu,U) + ex::div(mu*ex::faceDotGrad(U)));
        }

        io.write<colocated>();

        #include "colocatedContinuityErrors.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
}
