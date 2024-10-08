#include "Hysing.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "exSchemes.H"
#include <fstream>

#include "meshFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

List<scalar> time;
List<scalar> height;
List<scalar> velocity;

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

Hysing::~Hysing()
{}

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
    const colocatedScalarDirection& alpha =
        runTime_.lookupObjectRef<colocatedScalarField>("alpha")[0][0];

    const staggeredScalarField& Us =
        runTime_.lookupObjectRef<staggeredScalarField>("U");

    const colocatedVectorDirection U(ex::reconstruct(Us)()[0][0]);

    const colocatedVectorDirection& cc =
        alpha.fvMsh().metrics<colocated>().cellCenters()[0][0];

    const colocatedScalarDirection& cv =
        alpha.fvMsh().metrics<colocated>().cellVolumes()[0][0];

    const scalar v(gSum(cv*alpha));
    const vector h(gSum(cv*alpha*cc)/v);
    const vector u(gSum(cv*alpha*U)/v);

    if (Pstream::master())
        filePtr_()
            << runTime_.time().value()
            << " " << h.y() << " " << u.y() << nl;

    return true;
}

}

}

}

}
