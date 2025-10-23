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

// Runge-Kutta stage sources

FastPtrList<colocatedVectorField> stageSourcesA
(
    rk.stageList<colocatedVectorField>("A")
);

FastPtrList<colocatedVectorField> stageSourcesB
(
    rk.stageList<colocatedVectorField>("B")
);

linearSystem<stencil,vector,colocated> USysA("A(U)", U);
linearSystem<stencil,vector,colocated> USysB("B(U)", U);

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

colocatedScalarFaceField vf
(
    "vf",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

vf.setRestrictionScheme("areaWeighted");
