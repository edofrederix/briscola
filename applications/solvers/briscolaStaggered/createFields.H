// Velocity

staggeredScalarField U
(
    "U",
    fvMsh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE,
    true,
    true
);

U = Zero;
U.correctBoundaryConditions();

// Reconstructed colocated velocity

colocatedVectorField Uc
(
    "Uc",
    fvMsh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE,
    true,
    true
);

Uc = Zero;

// Reconstruct Uc at every time-step,
// needed for some function objects

Switch colocatedReconstruction
(
    solverDict.lookupOrDefault
    (
        "colocatedReconstruction",
        false
    )
);

// Pressure

colocatedScalarField p
(
    "p",
    fvMsh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE,
    true,
    true
);

p = Zero;
p.correctBoundaryConditions();

// Viscosity

const scalar nu = readScalar(solverDict.lookup("nu"));

// Source

staggeredScalarField exSource
(
    "source",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

staggeredScalarField imSourceCoeff
(
    "sourceCoeff",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

exSource = scalarList(solverDict.lookup("source"));
imSourceCoeff = scalarList(solverDict.lookup("sourceCoeff"));

// Immersed boundary

if (fvMsh.immersedBoundaryPresent())
{
    // Set source terms to zero inside IB
    forAllCells(exSource,l,d,i,j,k)
    {
        forAll(fvMsh.IB<staggered>(), ib)
        {
            if (fvMsh.IB<staggered>()[ib].mask()(l,d,i,j,k))
            {
                exSource(l,d,i,j,k) = Zero;
                imSourceCoeff(l,d,i,j,k) = Zero;
            }
        }
    }
}

Switch massSource = false;

forAll(U.IBC(), ib)
{
    const word IBname = fvMsh.IB<staggered>()[ib].name();

    if
    (
        U.IBC()[ib].Jac()
        && U.subDict("boundaryConditions").subDict(IBname)
            .lookupOrDefault<Switch>("massSource", true)
    )
    {
        massSource = true;
    }
}

// Solvers

linearSystem<stencil,scalar,staggered> USys(U);

autoPtr<solver<stencil,scalar,staggered>> USolve
(
    solver<stencil,scalar,staggered>::New
    (
        USys.x().name(),
        fvMsh
    )
);

autoPtr<PoissonSolver<symmStencil,scalar,colocated>> Poisson
(
    PoissonSolver<symmStencil,scalar,colocated>::New
    (
        p.name(),
        fvMsh
    )
);

staggeredLowerFaceScalarField phi("phi", fvMsh);

phi = ex::faceFlux(U);

// Adams-Bashforth source

staggeredScalarField DivU
(
    "DivU",
    fvMsh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE,
    true,
    false
);

DivU = ex::div(phi,U);