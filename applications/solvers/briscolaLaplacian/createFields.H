// Temperature

meshField<scalar,MESHTYPE> T
(
    "T",
    fvMsh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE,
    true,
    true
);

T = Zero;
T.correctBoundaryConditions();

// Diffusivity

meshField<hexScalar,MESHTYPE> lambda
(
    "lambda",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

lambda = readScalar(solverDict.lookup("lambda"));

// Source

meshField<scalar,MESHTYPE> source
(
    "source",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

source = readScalar(solverDict.lookup("source"));

// System and solver

linearSystem<stencil,scalar,MESHTYPE> TSys(T);

autoPtr<solver<stencil,scalar,colocated>> TSolve
(
    solver<stencil,scalar,colocated>::New
    (
        TSys.x().name(),
        fvMsh
    )
);