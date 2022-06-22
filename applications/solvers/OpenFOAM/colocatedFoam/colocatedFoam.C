#include "fvCFD.H"
#include "pimpleControl.H"

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"

    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    while (pimple.run(runTime))
    {
        #include "CourantNo.H"
        #include "setDeltaT.H"

        const dimensionedScalar deltaT(runTime.deltaT());

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
          + fvc::grad(p)
        );

        solve(UEqn);

        U += deltaT*fvc::grad(p);
        U.correctBoundaryConditions();

        phi = fvc::flux(U);

        fvScalarMatrix pEqn
        (
            fvm::laplacian(p) == fvc::div(phi)/deltaT
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        U -= deltaT*fvc::grad(p);
        U.correctBoundaryConditions();

        phi -= deltaT*fvc::snGrad(p)*mesh.magSf();

        #include "continuityErrs.H"

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    return 0;
}