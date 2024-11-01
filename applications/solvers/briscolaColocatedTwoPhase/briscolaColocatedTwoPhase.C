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

        Info << "Time = " << runTime.timeName() << endl;

        U.setOldTime();
        p.setOldTime();

        // Update the two-phase model and specific volumes

        twoPhase.correct();

        v = 1.0/rho;
        vf = ex::interp(v);

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
                USys -= C*exSource*v;
                USys -= rk.stageSum(stageSourcesA, stageSourcesB);

                if (rk.implicitStageA())
                {
                    USysA = -im::div(phi,U);
                    USys -= A*USysA;
                }

                if (rk.implicitStageB())
                {
                    USysB =
                        v
                      * (
                            im::laplacian(mu,U)
                          + ex::div(mu*ex::faceFlux(T(ex::grad(U))))
                          + im::source(imSourceCoeff,U)
                        );

                    USys -= B*USysB;
                }

                USys -= C*twoPhase.g();

                // Solve predictor

                USolve->solve(USys);

                // Pressure equation

                const colocatedLowerFaceScalarField phiStar
                (
                    ex::faceFlux(U)
                  + C*deltaT*twoPhase.surfaceTension()*vf
                );

                Poisson->solve(p, ex::div(phiStar)/(-C*deltaT), vf);

                // Rhie-Chow correction

                U -=
                    C*deltaT
                  * ex::reconstruct
                    (
                        Poisson->flux()/vf
                      - twoPhase.surfaceTension()
                    )*v;

                U.correctBoundaryConditions();

                if (rk.lastStage())
                    phi = phiStar - C*deltaT*Poisson->flux();
            }

            // Store Runge-Kutta sources

            if (!rk.lastStage())
            {
                stageSourcesA[stage-1] =
                    rk.solve() && rk.implicitStageA()
                  ? USysA.evaluate()
                  : -ex::div(phi,U);

                stageSourcesB[stage-1] =
                    rk.solve() && rk.implicitStageB()
                  ? USysB.evaluate()
                  : v
                  * (
                        ex::laplacian(mu,U)
                      + ex::div(mu*ex::faceFlux(T(ex::grad(U))))
                      + ex::source(imSourceCoeff,U)
                    );
            }
        }

        io.write<colocated>();

        #include "colocatedContinuityErrors.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
}
