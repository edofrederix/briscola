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
  * vector(twoPhaseDict.lookupOrDefault("source", vector::zero));

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

// Runge-Kutta stage sources

FastPtrList<colocatedVectorField>& stageSourcesA =
    rk.newStageList<colocatedVectorField>("A");

FastPtrList<colocatedVectorField>& stageSourcesB =
    rk.newStageList<colocatedVectorField>("B");

tmp<colocatedVectorSystem> tUSysA;
tmp<colocatedVectorSystem> tUSysB;

// Specific volumes

colocatedScalarField v
(
    "v",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

v = Zero;

colocatedScalarFaceField vf
(
    "vf",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

vf.setRestrictionScheme("areaWeighted");
v = Zero;

// Reference to the viscosity

const colocatedScalarFaceField& mu = twoPhase.mu();
