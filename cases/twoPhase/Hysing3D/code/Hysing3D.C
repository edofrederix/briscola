#include "Hysing3D.H"
#include "addToRunTimeSelectionTable.H"
#include "meshFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

defineTypeNameAndDebug(Hysing3D, 0);

addToRunTimeSelectionTable
(
    functionObject,
    Hysing3D,
    dictionary
);

bool Hysing3D::read(const dictionary& dict)
{
    if (runTime_.time().value() == 0.0)
    {
        colocatedScalarField& alpha =
            runTime_.lookupObjectRef<colocatedScalarField>("alpha");

        alpha = Zero;

        const colocatedVertexVectorField& v =
            alpha.fvMsh().metrics<colocated>().vertexCenters();

        forAllCells(alpha, i, j, k)
            for (int vi = 0; vi < 8; vi++)
                alpha(i,j,k) +=
                    0.125
                  * (
                        mag(v(i,j,k)[vi] - vector(0,0.5,0))
                     <= 0.25
                    );

        alpha.correctBoundaryConditions();
    }

    return true;
}

}

}

}

}
