#include "TaylorBubble.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
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

defineTypeNameAndDebug(TaylorBubble, 0);

addToRunTimeSelectionTable
(
    functionObject,
    TaylorBubble,
    dictionary
);

TaylorBubble::TaylorBubble
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    briscolaFunctionObject(name, runTime, dict)
{
    read(dict);
}

TaylorBubble::~TaylorBubble()
{}

bool TaylorBubble::read(const dictionary& dict)
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
                        mag(vector(v(i,j,k)[vi].x(), v(i,j,k)[vi].y(), 0))
                     <= 0.4
                    )
                  * (v(i,j,k)[vi].z() > 6.5 && v(i,j,k)[vi].z() < 9.5);

        alpha.correctBoundaryConditions();
    }

    return true;
}

bool TaylorBubble::execute()
{
    return true;
}

}

}

}

}
