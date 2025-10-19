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

staggeredScalarField imSourceCoeff
(
    "sourceCoeff",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

exSource =
    fvMsh.metrics<staggered>().fluidMask()
  * scalarList(solverDict.lookup("source"));

imSourceCoeff =
    fvMsh.metrics<staggered>().fluidMask()
  * scalarList(solverDict.lookup("sourceCoeff"));

// Runge-Kutta stage sources

FastPtrList<staggeredScalarField> stageSourcesA
(
    rk.stageList<staggeredScalarField>("A")
);

FastPtrList<staggeredScalarField> stageSourcesB
(
    rk.stageList<staggeredScalarField>("B")
);

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

staggeredScalarFaceField phi("phi", fvMsh);

phi = ex::faceFlux(U);
