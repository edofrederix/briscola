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

phi.setRestrictionScheme("flux");
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

colocatedVectorField exSource
(
    "source",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

colocatedScalarField imSourceCoeff
(
    "sourceCoeff",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

exSource = vector(solverDict.lookupOrDefault("source", vector::zero));
imSourceCoeff = scalar(solverDict.lookupOrDefault("sourceCoeff", 0.0));

// Solvers

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

// Stencils

linearSystem<stencil,vector,colocated> LapU(U);

// Adams-Bashforth source

colocatedVectorField DivU
(
    "DivU",
    fvMsh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE,
    true,
    false
);

DivU = ex::div(phi,U);
