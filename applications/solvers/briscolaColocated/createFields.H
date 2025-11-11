// Velocity

colocatedVectorField U
(
    "U",
    fvMsh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE,
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
    true
);

p = Zero;
p.correctBoundaryConditions();

// Flux

colocatedScalarFaceField phi
(
    "phi",
    fvMsh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE,
    true
);

phi = Zero;
phi.setRestrictionScheme("flux");

// Viscosity

const scalar nu = readScalar(solverDict.lookup("nu"));

// Source

colocatedVectorField exSource
(
    "source",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

exSource =
    fvMsh.metrics<colocated>().fluidMask()
  * solverDict.lookupOrDefault("source", vector::zero);

// Runge-Kutta stage sources

FastPtrList<colocatedVectorField>& stageSourcesA =
    rk.newStageList<colocatedVectorField>("A");

FastPtrList<colocatedVectorField>& stageSourcesB =
    rk.newStageList<colocatedVectorField>("B");

tmp<colocatedVectorSystem> tUSysA;
tmp<colocatedVectorSystem> tUSysB;

// Solvers

linearSystem<stencil,vector,colocated> USys("NS(U)", U, true);

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
        p.name(),
        fvMsh
    )
);

Poisson->enableFlux();
