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

colocatedLowerFaceScalarField phi
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

colocatedLowerFaceScalarField nu
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

linearSystem<symmStencil,vector,colocated> USys(U);
linearSystem<symmStencil,scalar,colocated> pSys(p);

autoPtr<solver<symmStencil,vector,colocated>> USolve
(
    solver<symmStencil,vector,colocated>::New
    (
        USys.x().name(),
        fvMsh
    )
);

autoPtr<PoissonSolver<symmStencil,scalar,colocated>> Poisson
(
    PoissonSolver<symmStencil,scalar,colocated>::New
    (
        pSys.x().name(),
        fvMsh
    )
);

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
