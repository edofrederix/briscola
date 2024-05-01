#include "initialCondition.H"
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

List<scalar> time;
List<scalar> height;
List<scalar> velocity;

defineTypeNameAndDebug(initialCondition, 0);

addToRunTimeSelectionTable
(
    functionObject,
    initialCondition,
    dictionary
);

initialCondition::initialCondition
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

initialCondition::~initialCondition()
{}

bool initialCondition::read(const dictionary& dict)
{
    if (runTime_.time().value() == 0.0)
    {
        colocatedScalarField& alpha =
            runTime_.lookupObjectRef<colocatedScalarField>("alpha");

        alpha = Zero;

        const meshField<vertexVector,colocated>& vertex =
            alpha.fvMsh().metrics<colocated>().vertexCenters();


        forAllCells(alpha, i, j, k)
        {
            scalar R = 0.25;

            scalar x0 = vertex(i,j,k).lba().x();
            scalar x1 = vertex(i,j,k).rba().x();
            scalar y0 = vertex(i,j,k).lba().y() - 0.5;
            scalar y1 = vertex(i,j,k).lta().y() - 0.5;

            if (x0 < -1e-12)
            {
                scalar aux = x0;
                x0 = -x1;
                x1 = -aux;
            }

            if (y0 < -1e-12)
            {
                scalar aux = y0;
                y0 = -y1;
                y1 = -aux;
            }

            scalar totalVol = (x1 - x0) * (y1 - y0);

            if (Foam::sqr(x1) + Foam::sqr(y1) < Foam::sqr(R))
            {
                alpha(i,j,k) = 1;
            }
            else if (Foam::sqr(x0) + Foam::sqr(y1) < Foam::sqr(R))
            {
                if (Foam::sqr(x1) + Foam::sqr(y0) < Foam::sqr(R))
                {
                    scalar aux = Foam::sqrt(Foam::sqr(R) - Foam::sqr(y1));

                    scalar vol = 0.5 *
                        (
                            Foam::sqr(R) * Foam::asin(x1 / R)
                          - Foam::sqr(R) * Foam::asin(aux / R)
                          + x1 * Foam::sqrt(Foam::sqr(R) - Foam::sqr(x1))
                          - aux * Foam::sqrt(Foam::sqr(R) - Foam::sqr(aux))
                        )
                      - y0 * (x1 - aux)
                      + (y1 - y0) * (aux - x0);

                    alpha(i,j,k) = vol / totalVol;

                }
                else
                {
                    scalar vol = 0.5 *
                        (
                            Foam::sqr(R) * Foam::asin(y1 / R)
                          - Foam::sqr(R) * Foam::asin(y0 / R)
                          + y1 * Foam::sqrt(Foam::sqr(R) - Foam::sqr(y1))
                          - y0 * Foam::sqrt(Foam::sqr(R) - Foam::sqr(y0))
                        )
                      - x0 * (y1 - y0);

                    alpha(i,j,k) = vol / totalVol;
                }
            }
            else if (Foam::sqr(x1) + Foam::sqr(y0) < Foam::sqr(R))
            {
                scalar vol = 0.5 *
                    (
                        Foam::sqr(R) * Foam::asin(x1 / R)
                      - Foam::sqr(R) * Foam::asin(x0 / R)
                      + x1 * Foam::sqrt(Foam::sqr(R) - Foam::sqr(x1))
                      - x0 * Foam::sqrt(Foam::sqr(R) - Foam::sqr(x0))
                    )
                  - y0 * (x1 - x0);

                alpha(i,j,k) = vol / totalVol;
            }
            else if (Foam::sqr(x0) + Foam::sqr(y0) < Foam::sqr(R))
            {
                scalar aux = Foam::sqrt(Foam::sqr(R) - Foam::sqr(y0));

                scalar vol = 0.5 *
                    (
                        Foam::sqr(R) * Foam::asin(aux / R)
                      - Foam::sqr(R) * Foam::asin(x0 / R)
                      + aux * Foam::sqrt(Foam::sqr(R) - Foam::sqr(aux))
                      - x0 * Foam::sqrt(Foam::sqr(R) - Foam::sqr(x0))
                    )
                  - y0 * (aux - x0);

                alpha(i,j,k) = vol / totalVol;

            }
            else
            {
                alpha(i,j,k) = 0;
            }
        }

        alpha.correctBoundaryConditions();

        time.clear();
        height.clear();
        velocity.clear();

        time.append(runTime_.time().value());
        height.append(0.5);
        velocity.append(0.0);
    }

    return true;
}

bool initialCondition::execute()
{
    const colocatedScalarField& alpha =
        runTime_.lookupObjectRef<colocatedScalarField>("alpha");

    const colocatedVectorField& U =
        runTime_.lookupObjectRef<colocatedVectorField>("U");

    const meshField<vector,colocated>& cc =
        alpha.fvMsh().metrics<colocated>().cellCenters();

    const meshField<scalar,colocated>& cv =
        alpha.fvMsh().metrics<colocated>().cellVolumes();

    scalar vol = 0;
    scalar pos = 0;
    scalar vel = 0;

    forAllCells(alpha, i, j, k)
    {
        vol += cv(i,j,k) * alpha(i,j,k);
        pos += cv(i,j,k) * alpha(i,j,k) * cc(i,j,k)[1];
        vel += cv(i,j,k) * alpha(i,j,k) * U(i,j,k)[1];
    }

    reduce(vol, sumOp<scalar>());
    reduce(pos, sumOp<scalar>());
    reduce(vel, sumOp<scalar>());

    pos /= vol;
    vel /= vol;

    time.append(runTime_.time().value());
    height.append(pos);
    velocity.append(vel);

    return true;
}

bool initialCondition::end()
{
    const colocatedScalarField& alpha =
        runTime_.lookupObjectRef<colocatedScalarField>("alpha");

    const colocatedVectorField& U =
        runTime_.lookupObjectRef<colocatedVectorField>("U");

    const meshField<vector,colocated>& cc =
        alpha.fvMsh().metrics<colocated>().cellCenters();

    const meshField<scalar,colocated>& cv =
        alpha.fvMsh().metrics<colocated>().cellVolumes();

    scalar vol = 0;
    scalar pos = 0;
    scalar vel = 0;

    forAllCells(alpha, i, j, k)
    {
        vol += cv(i,j,k) * alpha(i,j,k);
        pos += cv(i,j,k) * alpha(i,j,k) * cc(i,j,k)[1];
        vel += cv(i,j,k) * alpha(i,j,k) * U(i,j,k)[1];
    }

    reduce(vol, sumOp<scalar>());
    reduce(pos, sumOp<scalar>());
    reduce(vel, sumOp<scalar>());

    pos /= vol;
    vel /= vol;

    time.append(runTime_.time().value());
    height.append(pos);
    velocity.append(vel);

    if (Pstream::myProcNo() == 0)
    {
        std::ofstream outfile;
        outfile.open("hysingResults.txt", std::ios_base::app);

        int n = time.size();

        for (int i = 0; i < n; i++)
        {
            outfile << time[i] << " " << height[i] << " " << velocity[i] << "\n";
        }
    }

    return true;
}

}

}

}

}
