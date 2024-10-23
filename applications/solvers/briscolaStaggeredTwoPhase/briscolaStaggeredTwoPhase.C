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
        const scalar deltaT0 = runTime.deltaT0Value();

        Info << "Time = " << runTime.timeName() << endl;

        U.setOldTime();
        p.setOldTime();

        // Update the two-phase model and specific volumes

        twoPhase.correct();

        v = 1.0/rho;
        vcf = ex::coloFaceInterp(v);

        // Predictor

        USys = im::ddt(U);

        USys -= im::source(imSourceCoeff,U);
        USys -= v*exSource;

        USys -= im::laplacian(mu,U,0.5)*v;
        USys -= ex::laplacian(mu,U,0.5)*v;

        USys -= 0.5*(deltaT/deltaT0)*H;

        phi = ex::faceFlux(U);

        H = ex::div(phi,U);

        USys += (1.0 + 0.5*(deltaT/deltaT0))*H;
        USys -= twoPhase.gravity();
        USys -= (ex::grad(mu) & stagT(ex::grad(U)))*v;

        // Solve predictor

        USolve->solve(USys);

        U += deltaT*ex::stagReconstruct(twoPhase.surfaceTension())*v;
        U.correctBoundaryConditions();

        // Pressure equation

        Poisson->solve(p, ex::coloDiv(U)/(-deltaT), vcf);

        // Correct velocity

        U -= deltaT*ex::stagReconstruct(Poisson->flux()/vcf)*v;
        U.correctBoundaryConditions();

        io.write<colocated>();
        io.write<staggered>();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
}
