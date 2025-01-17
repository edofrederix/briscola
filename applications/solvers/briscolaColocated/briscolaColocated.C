#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"

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
            "briscolaColocatedDict",
            fvMsh.time().system(),
            fvMsh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

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

                if (rk.imStageA())
                {
                    USysA = -im::div(phi,U);
                    USys -= A*USysA;
                }

                if (rk.imStageB())
                {
                    USysB =
                        im::laplacian(nu,U)
                      + im::source(imSourceCoeff,U);

                    USys -= B*USysB;
                }

                // Solve predictor

                USolve->solve(USys);

                // Pressure equation

                const colocatedFaceScalarField phiStar(ex::faceFlux(U));

                Poisson->solve(p, ex::div(phiStar)/(-C*deltaT));

                // Rhie-Chow correction

                U -= C*deltaT*ex::grad(p);
                U.correctBoundaryConditions();

                if (rk.lastStage())
                    phi = phiStar - C*deltaT*Poisson->flux();
            }

            // Store Runge-Kutta sources

            if (!rk.lastStage())
            {
                stageSourcesA[stage-1] =
                    rk.solve() && rk.imStageA()
                  ? USysA.evaluate()
                  : -ex::div(phi,U);

                stageSourcesB[stage-1] =
                    rk.solve() && rk.imStageB()
                  ? USysB.evaluate()
                  : ex::laplacian(nu,U) + ex::source(imSourceCoeff,U);
            }
        }

        io.write<colocated>();

        #include "colocatedContinuityErrors.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
}
