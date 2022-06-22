#include "argList.H"
#include "timeSelector.H"
#include "volFields.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    const vectorField& C = mesh.C();
    vectorField& u = U.field();

    forAll(u, celli)
    {
        u[celli].x() =
            5300.0/360.0*2.0
          * (
                (1.0-Foam::sqr(C[celli].y()))
              + 0.4*Foam::sin(C[celli].z()*3)
              + 0.1*Foam::sin(C[celli].y()*3.1415927)
            );

        u[celli].z() =
            5300.0/360.0*2.0
          * Foam::sin(C[celli].x()*3)
          * Foam::sin(C[celli].y()*3.1415927);
    }

    U.correctBoundaryConditions();
    U.write();

    return 0;
}


// ************************************************************************* //
