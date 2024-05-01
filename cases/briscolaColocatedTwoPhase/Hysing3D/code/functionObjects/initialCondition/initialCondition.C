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
            for (int aux = 0; aux < 8; aux++)
            {
                alpha(i,j,k) += 0.125 * (Foam::mag(vertex(i,j,k)[aux] - vector(0,0.5,0)) < 0.25);
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
