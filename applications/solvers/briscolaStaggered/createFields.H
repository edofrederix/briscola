// Velocity

staggeredScalarField U
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

// Reconstructed colocated velocity

colocatedVectorField Uc
(
    "Uc",
    fvMsh,
    IOobject::NO_READ,
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
    true,
    true
);

p = Zero;
p.correctBoundaryConditions();

// Viscosity

staggeredFaceScalarField nu
(
    "nu",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

nu = readScalar(solverDict.lookup("nu"));

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

exSource = scalarList(solverDict.lookup("source"));
imSourceCoeff = scalarList(solverDict.lookup("sourceCoeff"));

// Immersed boundary

immersedBoundary<scalar,staggered> IBs(solverDict, fvMsh);

if (solverDict.found("ImmersedBoundary"))
{
    forAllLevels(exSource,l,d,i,j,k)
    {
        if (IBs.mask()(l,d,i,j,k) == 1)
        {
            // Set source terms to zero inside IB
            exSource = Zero;
            imSourceCoeff = Zero;
        }
    }
}

// Solvers

linearSystem<stencil,scalar,staggered> USys(U);

autoPtr<solver<stencil,scalar,staggered>> USolve
(
    solver<stencil,scalar,staggered>::New
    (
        USys.x().name(),
        fvMsh,
        (
            (IBs.type() == "penalization" || IBs.type() == "Fadlun") ?
            nullptr : &IBs
        )
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

staggeredFaceScalarField phi("phi", fvMsh);

phi = ex::faceFlux(U);

// Stencils

linearSystem<stencil,scalar,staggered> LapU(U);

// Adams-Bashforth source

staggeredScalarField DivU
(
    "DivU",
    fvMsh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE,
    true,
    false
);

DivU = ex::div(phi,U);
