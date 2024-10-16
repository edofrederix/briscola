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

// Pressure gradient

colocatedVectorField G
(
    "grad(p)",
    fvMsh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE,
    true
);

G = Zero;

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

exSource = vector(twoPhaseDict.lookupOrDefault("source", vector::zero));
imSourceCoeff = scalar(twoPhaseDict.lookupOrDefault("sourceCoeff", 0.0));

// Solvers

linearSystem<stencil,vector,colocated> USys(U);
linearSystem<symmStencil,scalar,colocated> pSys(p);

autoPtr<solver<stencil,vector,colocated>> USolve
(
    solver<stencil,vector,colocated>::New
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

colocatedLowerFaceScalarField vf
(
    "vf",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

vf.setRestrictionScheme("faceAreaWeighted");

// Adams-Bashforth source

colocatedVectorField H
(
    "H",
    fvMsh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE,
    true
);

H = Zero;
