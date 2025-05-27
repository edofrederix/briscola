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

                if (rk.imStageA())
                {
                    USysA = -im::div(phi,U);
                    USys -= A*USysA;
                }

                if (rk.imStageB())
                {
                    USysB =
                        v
                      * (
                            im::laplacian(mu,U)
                          + (ex::grad(mu) & stagT(ex::grad(U)))
                          + im::source(imSourceCoeff,U)
                        );

                    USys -= B*USysB;
                }

                USys -= C*twoPhase.buoyancy()*v;

                // Solve predictor

                USolve->solve(USys);

                U += C*deltaT*ex::stagReconstruct(twoPhase.flux())*v;
                U.correctBoundaryConditions();

                // Pressure equation

                Poisson->solve(p, ex::coloDiv(U)/(-C*deltaT), vcf);

                // Correct velocity

                U -= C*deltaT*ex::stagReconstruct(Poisson->flux()/vcf)*v;
                U.correctBoundaryConditions();
            }

            // Store Runge-Kutta sources

            if (rk.storeStageA())
                stageSourcesA[stage-1] =
                    rk.solve() && rk.imStageA()
                  ? USysA.evaluate()
                  : -ex::div(phi,U);

            if (rk.storeStageB())
                stageSourcesB[stage-1] =
                    rk.solve() && rk.imStageB()
                  ? USysB.evaluate()
                  : v
                  * (
                        ex::laplacian(mu,U)
                      + (ex::grad(mu) & stagT(ex::grad(U)))
                      + ex::source(imSourceCoeff,U)
                    );
        }

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
