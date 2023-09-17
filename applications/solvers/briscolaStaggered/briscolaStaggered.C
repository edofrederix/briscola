#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"

#include "fv.H"
#include "immersedBoundary.H"

using namespace Foam;
using namespace briscola;
using namespace fv;
using namespace ibm;

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"
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

        // Predictor

        USys = im::ddt(U);

        USys -= im::source(imSourceCoeff,U);
        USys -= exSource;

        LapU = im::laplacian(nu,U);
        USys -= 0.5*LapU;
        USys -= 0.5*LapU.evaluate();

        USys -= 0.5*DivU;
        phi = ex::faceFlux(U);
        DivU = ex::div(phi,U);
        USys += 1.5*DivU;

        USys.correctBoundaries();

        // Immersed boundary

        if (solverDict.found("ImmersedBoundary"))
        {
            IBs.penalization(USys);

            if
            (
                solverDict.subDict("ImmersedBoundary")
                    .lookupOrDefault("IBM", true)
            )
            {
                IBs.IBM(USys);
                Info << "IBM!" << endl;
            }
        }

        // Solve predictor

        USolve->solve(USys);

        // Pressure equation

        if
        (
            solverDict.found("ImmersedBoundary")
            && solverDict.subDict("ImmersedBoundary")
                .lookupOrDefault("PCorr", false)
        )
        {
            for (int iter = 0; iter < 3; iter++)
            {
                Poisson->solve(p, -ex::coloDiv(U)/deltaT + IBc.exPCorr(p));
            }
            Info << "PCorr!" << endl;
        }
        else
        {
            Poisson->solve(p, -ex::coloDiv(U)/deltaT);
        }

        // Correction

        U -= deltaT*ex::stagGrad(p);
        U.correctBoundaryConditions();

        if (fvMsh.time().writeTime())
            Uc = ex::reconstruct(U);

        io.write<colocated>();
        io.write<staggered>();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
}
