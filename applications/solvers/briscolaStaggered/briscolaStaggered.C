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

        // Predictor

        USys = im::ddt(U);

        USys -= im::source(imSourceCoeff,U);
        USys -= exSource;

        USys -= im::laplacian(0.5*nu,U);
        USys -= ex::laplacian(0.5*nu,U);

        USys -= 0.5*H;
        phi = ex::faceFlux(U);
        H = ex::div(phi,U);
        USys += 1.5*H;

        // Solve predictor

        USolve->solve(USys);

        // Pressure equation

        Poisson->solve(p, IBMCorr(ex::coloDiv(U),U)/(-deltaT));

        // Correction

        U -= deltaT*ex::stagGrad(p);
        U.correctBoundaryConditions();

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
