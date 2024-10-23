// Temperature

meshField<scalar,colocated> T
(
    "T",
    fvMsh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE,
    true
);

T = Zero;
T.correctBoundaryConditions();
T.correctImmersedBoundaryConditions();

// Diffusivity

meshField<lowerFaceScalar,colocated> lambda
(
    "lambda",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

lambda = readScalar(solverDict.lookup("lambda"));

// Source

meshField<scalar,colocated> source
(
    "source",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

source = readScalar(solverDict.lookup("source"));

// System and solver

autoPtr<PoissonSolver<symmStencil,scalar,colocated>> Poisson
(
    PoissonSolver<symmStencil,scalar,colocated>::New
    (
        T.name(),
        fvMsh
    )
);
