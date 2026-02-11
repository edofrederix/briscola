#include "bubblePipe2D.H"
#include "addToRunTimeSelectionTable.H"
#include "List.H"
#include "PstreamReduceOps.H"
#include <fstream>
#include "Time.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

defineTypeNameAndDebug(bubblePipe2D, 0);

addToRunTimeSelectionTable
(
    functionObject,
    bubblePipe2D,
    dictionary
);

bubblePipe2D::bubblePipe2D
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    briscolaFunctionObject(name, runTime, dict),
    fvMsh_(runTime.lookupObject<fvMesh>("briscolaMeshDict"))
{
    read(dict);
}

bubblePipe2D::~bubblePipe2D()
{}

bool bubblePipe2D::read(const dictionary& dict)
{
    const vector u(0, 0.2, 0);

    if (runTime_.time().value() == 0.0)
    {
        if (runTime_.foundObject<staggeredScalarField>("U"))
        {
            staggeredScalarField& U =
                runTime_.lookupObjectRef<staggeredScalarField>("U");

            U = Zero;

            const tensor base =
                U.fvMsh().msh().cast<rectilinearMesh>().base();

            forAllCells(U,d,i,j,k)
                U(d,i,j,k) = (base & u)[d];

            U.correctBoundaryConditions();
        }
        else
        {
            colocatedVectorField& U =
                runTime_.lookupObjectRef<colocatedVectorField>("U");

            U = Zero;

            forAllCells(U,i,j,k)
                U(i,j,k) = u;

            U.correctBoundaryConditions();
        }

        colocatedScalarField& alpha =
            runTime_.lookupObjectRef<colocatedScalarField>("alpha");

        alpha = Zero;

        const colocatedVertexVectorField& v =
            fvMsh_.metrics<colocated>().vertexCenters();

        vectorList bubbleCenters(0);

        const scalar bubbleRadius(0.0025);
        const scalar pipeRadius(0.02);
        const scalar pipeLength(0.1);

        const int nBubbles(30);

        // std::srand(std::time({}));
        std::srand(123456789);

        if (Pstream::master())
        {
            while (bubbleCenters.size() < nBubbles)
            {
                scalar x =
                    (2.0*std::rand()/RAND_MAX - 1.0)
                  * (pipeRadius-bubbleRadius);

                scalar y = 0.0;
                scalar z =
                    bubbleRadius
                  + std::rand()/(double)RAND_MAX
                  * (pipeLength - 2.0 * bubbleRadius);

                vector newBubble(x,y,z);

                bool overlap = false;

                // Check overlap with other bubbles
                forAll(bubbleCenters, i)
                {
                    scalar distance = mag(bubbleCenters[i] - newBubble);

                    if (distance < 2.0*bubbleRadius + 1e-3)
                    {
                        overlap = true;
                        break;
                    }
                }

                if (!overlap)
                    bubbleCenters.append(newBubble);
            }
        }

        Pstream::scatter(bubbleCenters);

        forAllCells(alpha, i, j, k)
            for (int vi = 0; vi < 8; vi++)
                forAll(bubbleCenters, b)
                    alpha(i,j,k) +=
                        0.125
                      * (
                            mag
                            (
                                vector(v(i,j,k)[vi].x(), 0.0, v(i,j,k)[vi].z())
                              - bubbleCenters[b]
                            )
                         <= bubbleRadius
                        );

        alpha.correctBoundaryConditions();
    }

    return true;
}

}

}

}

}
