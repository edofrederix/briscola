#include "TaylorBubbleJSI.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

defineTypeNameAndDebug(TaylorBubbleJSI, 0);

addToRunTimeSelectionTable
(
    functionObject,
    TaylorBubbleJSI,
    dictionary
);

bool TaylorBubbleJSI::read(const dictionary& dict)
{
    const scalar bubbleHeight = 0.03;
    const scalar bubbleCenter = 0.124;
    const scalar bubbleRadius = 0.005;

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
                     <= bubbleRadius
                    )
                  * (
                        v(i,j,k)[vi].z() > (bubbleCenter - bubbleHeight/2)
                     && v(i,j,k)[vi].z() < (bubbleCenter + bubbleHeight/2)
                    );

        alpha.correctBoundaryConditions();
    }

    return true;
}

bool TaylorBubbleJSI::execute()
{
    return true;
}

}

}

}

}
