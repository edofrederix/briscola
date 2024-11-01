#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"
#include "immersedBoundaryConditionStaggeredMassSource.C"

#include "fv.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"
    #include "createRungeKuttaScheme.H"
    #include "createTimeControls.H"

    // Solver dictionary

    IOdictionary solverDict
    (
        IOobject
        (
            "briscolaStaggeredDict",
            fvMsh.time().system(),
            fvMsh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

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
                USys -= C*exSource;
                USys -= rk.stageSum(stageSourcesA, stageSourcesB);

                if (rk.implicitStageA())
                {
                    USysA = -im::div(phi,U);
                    USys -= A*USysA;
                }

                if (rk.implicitStageB())
                {
                    USysB =
                        im::laplacian(nu,U)
                      + im::source(imSourceCoeff,U);

                    USys -= B*USysB;
                }

                // Solve predictor

                USolve->solve(USys);

                // Pressure equation

                Poisson->solve(p, ibmCorr(ex::coloDiv(U),U)/(-C*deltaT));

                // Correction

                U -= C*deltaT*ex::stagReconstruct(Poisson->flux());
                U.correctBoundaryConditions();
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
                  : ex::laplacian(nu,U) + ex::source(imSourceCoeff,U);
            }
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
