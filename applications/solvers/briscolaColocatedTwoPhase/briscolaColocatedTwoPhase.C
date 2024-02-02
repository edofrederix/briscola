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

    Switch split = args.optionFound("split");

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

        v = icoTwoPhase.v<colocated>();
        vf = ex::interp(v);

        // Predictor, Eq. (A.1) of Dodd & Ferrante (2014)

        USys = im::ddt(U);

        USys -= im::laplacian(mu,U,0.5)*v;
        USys -= ex::laplacian(mu,U,0.5)*v;

        USys -= 0.5*(deltaT/deltaT0)*H;

        colocatedVectorField surfTen
        (
            ex::reconstruct(fa*icoTwoPhase.surfaceTension().potential())
        );

        H = ex::div(phi,U) - (ex::grad(mu) & ex::grad(U))*v;

        USys += (1.0 + 0.5*(deltaT/deltaT0))*H;

        USys -= icoTwoPhase.g();
        USys += ex::grad(p)*v - surfTen;

        // Solve predictor

        USolve->solve(USys);

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

            Poisson->solve(p, minRho*ex::div(phi)/(-deltaT) - ex::div(corr));

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

            phi -= deltaT*(ex::faceGrad(p) - ex::faceGrad(p.oldTime()))*fa*vf;
        }

        io.write<colocated>();

        #include "colocatedContinuityErrors.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
}
