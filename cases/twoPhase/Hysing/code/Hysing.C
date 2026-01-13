#include "Hysing.H"
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

defineTypeNameAndDebug(Hysing, 0);

addToRunTimeSelectionTable
(
    functionObject,
    Hysing,
    dictionary
);

Hysing::Hysing
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    briscolaFunctionObject(name, runTime, dict)
{
    read(dict);

    if (Pstream::master())
    {
        filePtr_.reset(new OFstream("data.txt"));
    }
}

bool Hysing::read(const dictionary& dict)
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
                        mag(vector(v(i,j,k)[vi].x(), v(i,j,k)[vi].y()-0.5, 0))
                     <= 0.25
                    );

        alpha.correctBoundaryConditions();
    }

    return true;
}

bool Hysing::execute()
{
    const colocatedScalarField& alpha =
        runTime_.lookupObject<colocatedScalarField>("alpha");

    const colocatedVectorField& U =
        runTime_.foundObject<colocatedVectorField>("U")
      ? runTime_.lookupObject<colocatedVectorField>("U")
      : runTime_.lookupObject<colocatedVectorField>("Uc");

    const colocatedVectorField& cc =
        alpha.fvMsh().metrics<colocated>().cellCenters();

    const colocatedScalarField& cv =
        alpha.fvMsh().metrics<colocated>().cellVolumes();

    const scalar v(gSum(cv*alpha)[0]);
    const vector h(gSum(cv*alpha*cc)[0]/v);
    const vector u(gSum(cv*alpha*U)[0]/v);

    if (Pstream::master())
        filePtr_()
            << runTime_.time().value()
            << " " << h.y() << " " << u.y() << " " << v << nl;

    return true;
}

}

}

}

}
