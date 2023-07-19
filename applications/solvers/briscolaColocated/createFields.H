// Velocity

colocatedVectorField U
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

// Flux

colocatedFaceScalarField phi
(
    "phi",
    fvMsh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE,
    true
);

phi = ex::faceFlux(U);

// Viscosity

colocatedFaceScalarField nu
(
    "nu",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

nu = readScalar(solverDict.lookup("nu"));

// Source

colocatedVectorField source
(
    "source",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

colocatedScalarField sourceCoeff
(
    "sourceCoeff",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

source = vector(solverDict.lookupOrDefault("source", vector::zero));
sourceCoeff = scalar(solverDict.lookupOrDefault("sourceCoeff", 0.0));

linearSystem<stencil,vector,colocated> USys(U);
linearSystem<stencil,scalar,colocated> pSys(p);

autoPtr<solver<stencil,vector,colocated>> USolve
(
    solver<stencil,vector,colocated>::New
    (
        USys.x().name(),
        fvMsh
    )
);

autoPtr<PoissonSolver<stencil,scalar,colocated>> Poisson
(
    PoissonSolver<stencil,scalar,colocated>::New
    (
        pSys.x().name(),
        fvMsh
    )
);
