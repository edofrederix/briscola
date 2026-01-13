#include "slottedDisk.H"
#include "addToRunTimeSelectionTable.H"
#include "meshFields.H"
#include "uniformMesh.H"
#include "exSchemesFaceFlux.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

defineTypeNameAndDebug(slottedDisk, 0);

addToRunTimeSelectionTable
(
    functionObject,
    slottedDisk,
    dictionary
);

tmp<colocatedScalarField> slottedDisk::alphaRef() const
{
    const label Nx = N_.x();
    const label Ny = N_.y();
    const label Nz = N_.z();

    const fvMesh& fvMsh =
        runTime_.lookupObject<fvMesh>("briscolaMeshDict");

    const meshField<vertexVector,colocated>& vertex =
        fvMsh.metrics<colocated>().vertexCenters();

    const uniformMesh& msh = fvMsh.msh().cast<uniformMesh>();

    const tensor base = msh.base();
    const vector s = (base.T() & msh.cellSize());

    // Reference discretization cell

    vertexVector ref;

    for (label k = 0; k < 2; k++)
        for (label j = 0; j < 2; j++)
            for (label i = 0; i < 2; i++)
                ref[k*4+j*2+i] =
                    vector(i*s.x()/Nx, j*s.y()/Ny, k*s.z()/Nz);

    // Compute the reference alpha field using a simple cell discretization

    tmp<colocatedScalarField> tAlpha =
        colocatedScalarField::New("alpha", fvMsh);

    colocatedScalarField& alpha = tAlpha.ref();

    alpha = Zero;

    forAllCells(alpha, i, j, k)
    {
        const labelVector ijk(i,j,k);

        for (label ii = 0; ii < Nx; ii++)
        for (label jj = 0; jj < Ny; jj++)
        for (label kk = 0; kk < Nz; kk++)
        {
            const vertexVector vert = vertex(ijk);
            const vector start = min(vert.lba(), vert.rtf());

            const vertexVector cell =
                ref
              + start
              + vector(ii*s.x()/Nx, jj*s.y()/Ny, kk*s.z()/Nz);

            for (label vi = 0; vi < 8; vi++)
                alpha(ijk) +=
                    0.125/(Nx*Ny*Nz)*inside(cell[vi]);
        }
    }

    return tAlpha;
}

slottedDisk::slottedDisk
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    briscolaFunctionObject(name, runTime, dict),
    T_(1.0),
    R_(0.15),
    C_(0.5,0.75,0.0),
    N_(8,8,1)
{
    read(dict);
}

bool slottedDisk::read(const dictionary& dict)
{
    const scalar t = runTime_.time().value();
    const scalar dt = runTime_.deltaT().value();

    if (t == 0.0)
    {
        // Set velocity

        colocatedVectorField& U =
            runTime_.lookupObjectRef<colocatedVectorField>("U");

        colocatedScalarFaceField& phi =
            runTime_.lookupObjectRef<colocatedScalarFaceField>("phi");

        const fvMesh& fvMsh = U.fvMsh();

        const colocatedVectorField& cc =
            fvMsh.metrics<colocated>().cellCenters();

        forAllCells(U, i, j, k)
            U(i,j,k) = v(cc(i,j,k), t, dt);

        U.correctBoundaryConditions();
        phi = ex::faceFlux(U);

        // Initialize volume fraction

        colocatedScalarField& alpha =
            runTime_.lookupObjectRef<colocatedScalarField>("alpha");

        alpha = alphaRef();
        alpha.correctBoundaryConditions();
    }

    return true;
}

bool slottedDisk::execute()
{
    const scalar t = runTime_.time().value();
    const scalar dt = runTime_.deltaT().value();

    // Update velocity

    colocatedVectorField& U =
        runTime_.lookupObjectRef<colocatedVectorField>("U");

    colocatedScalarFaceField& phi =
        runTime_.lookupObjectRef<colocatedScalarFaceField>("phi");

    const fvMesh& fvMsh = U.fvMsh();

    const colocatedVectorField& cc =
        fvMsh.metrics<colocated>().cellCenters();

    forAllCells(U, i, j, k)
        U(i,j,k) = v(cc(i,j,k), t, dt);

    U.correctBoundaryConditions();
    phi = ex::faceFlux(U);

    return true;
}

bool slottedDisk::end()
{
    const scalar t = runTime_.time().value();

    if (mag(t - T_) < 1e-12)
    {
        colocatedScalarField& alpha =
            runTime_.lookupObjectRef<colocatedScalarField>("alpha");

        const fvMesh& fvMsh = alpha.fvMsh();

        const colocatedScalarField& cv =
            fvMsh.metrics<colocated>().cellVolumes();

        const colocatedScalarField ref(alphaRef());

        scalar error = 0.0;
        scalar vol = 0.0;
        scalar volRef = 0.0;

        forAllCells(alpha, i, j, k)
        {
            labelVector ijk(i,j,k);

            error += cv(ijk)*mag(ref(ijk) - alpha(ijk));
            vol += cv(ijk)*alpha(ijk);
            volRef += cv(ijk)*ref(ijk);
        }

        reduce(error, sumOp<scalar>());
        reduce(vol, sumOp<scalar>());
        reduce(volRef, sumOp<scalar>());

        Info<< "Reference volume = " << volRef << nl
            << "Final volume = " << vol << nl
            << "Volume error = " << (volRef - vol)/vol << nl
            << "Shape error = " << error/vol << endl;
    }

    return true;
}

}

}

}

}
