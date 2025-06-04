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

exSource =
    fvMsh.metrics<colocated>().fluidMask()
  * vector(twoPhaseDict.lookupOrDefault("source", vector::zero));

imSourceCoeff =
    fvMsh.metrics<colocated>().fluidMask()
  * scalar(twoPhaseDict.lookupOrDefault("sourceCoeff", 0.0));

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

// Specific volumes

colocatedScalarField v
(
    "v",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

v.setRestrictionScheme("harmonicVolumeWeighted");

colocatedFaceScalarField vf
(
    "vf",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

vf.setRestrictionScheme("faceAreaWeighted");

// CNAB sources

colocatedVectorField H
(
    "H",
    fvMsh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE,
    true
);

H = ex::div(phi,U);

colocatedVectorField G
(
    "G",
    fvMsh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE,
    true
);

G = ex::laplacian(mu,U)*v;

// Modified or classic CNAB

const Switch modified(twoPhaseDict.lookupOrDefault<Switch>("modified", true));
