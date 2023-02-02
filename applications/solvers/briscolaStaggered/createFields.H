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
    IOobject::NO_READ,
    IOobject::AUTO_WRITE,
    true
);

Uc = Zero;

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

staggeredFaceScalarField nu
(
    "nu",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

nu = readScalar(solverDict.lookup("nu"));

// Source

staggeredScalarField source
(
    "source",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

source = scalarList(solverDict.lookup("source"));

linearSystem<stencil,scalar,staggered> USys(U);

autoPtr<solver<stencil,scalar,staggered>> USolve
(
    solver<stencil,scalar,staggered>::New
    (
        USys.x().name(),
        fvMsh
    )
);

autoPtr<PoissonSolver<stencil,scalar,colocated>> Poisson
(
    PoissonSolver<stencil,scalar,colocated>::New
    (
        p.name(),
        fvMsh
    )
);

staggeredFaceScalarField phi("phi", fvMsh);

phi = Zero;
