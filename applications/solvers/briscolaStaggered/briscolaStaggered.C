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
    #include "createBriscolaStaggeredMesh.H"
    #include "createRungeKuttaScheme.H"
    #include "createTimeControls.H"

    // Solver dictionary

    IOdictionary solverDict
    (
        IOobject
        (
            "briscolaSinglePhaseDict",
            fvMsh.time().system(),
            fvMsh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    #include "createFields.H"
    #include "createBriscolaIO.H"

    if (rk.imA())
        restrict(U);

    U.correctBoundaryConditions();

    while (runTime.run())
    {
        #include "staggeredCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        const scalar deltaT = runTime.deltaTValue();

        Info << "Time = " << runTime.name() << endl;

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

                if (rk.imStageA())
                {
                    tUSysA = -im::div(phi,U);
                    USys -= A*tUSysA.ref();
                }

                if (rk.imStageB())
                {
                    tUSysB = im::laplacian(nu,U);
                    USys -= B*tUSysB.ref();
                }

                // Solve predictor

                USolve->solve(USys);

                // Pressure equation

                Poisson->solve(p, ibmCorr(ex::coloDiv(U),U)/(-C*deltaT));

                // Correction

                U -= C*deltaT*ex::stagReconstruct(Poisson->flux());

                if (rk.imA())
                    restrict(U);

                U.correctBoundaryConditions();
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
                  : ex::laplacian(nu,U);
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
