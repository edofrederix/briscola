// Velocity

staggeredScalarField U
(
    "U",
    fvMsh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE,
    true
);

U = Zero;
U.correctBoundaryConditions();

// Reconstructed colocated velocity

colocatedVectorField Uc
(
    "Uc",
    fvMsh,
    IOobject::READ_IF_PRESENT,
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
    true
);

p = Zero;
p.correctBoundaryConditions();

// Flux

staggeredScalarFaceField phi
(
    "phi",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

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

// Project the source vector onto the local coordinate system

const tensor base =
    fvMsh.msh().cast<rectilinearMesh>().base();

exSource =
    fvMsh.metrics<staggered>().fluidMask()
  * list(base & solverDict.lookupOrDefault("source", vector::zero));

// Runge-Kutta stage sources

FastPtrList<staggeredScalarField>& stageSourcesA =
    rk.newStageList<staggeredScalarField>("A");

FastPtrList<staggeredScalarField>& stageSourcesB =
    rk.newStageList<staggeredScalarField>("B");

linearSystem<stencil,scalar,staggered> USysA("A(U)", U);
linearSystem<stencil,scalar,staggered> USysB("B(U)", U);

// Solvers

linearSystem<stencil,scalar,staggered> USys("NS(U)", U, true);

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

Poisson->enableFlux();
