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

exSource = scalarList(solverDict.lookup("source"));
imSourceCoeff = scalarList(solverDict.lookup("sourceCoeff"));

// Immersed boundary

if (fvMsh.ibs<staggered>().size())
{
    // Set source terms to zero inside IB
    forAll(fvMsh.ibs<staggered>(), ib)
    {
        forAllCells(exSource,l,d,i,j,k)
        {
            if (fvMsh.ibs<staggered>()[ib].mask()(l,d,i,j,k))
            {
                exSource(l,d,i,j,k) = Zero;
                imSourceCoeff(l,d,i,j,k) = Zero;
            }
        }
    }
}

// Solvers

linearSystem<stencil,scalar,staggered> USys(U, true);

autoPtr<solver<stencil,scalar,staggered>> USolve
(
    solver<stencil,scalar,staggered>::New
    (
        USys.x().name(),
        fvMsh
    )
);

autoPtr<PoissonSolver<symmStencil,scalar,colocated>> Poisson
(
    PoissonSolver<symmStencil,scalar,colocated>::New
    (
        p.name(),
        fvMsh
    )
);

Poisson->enableFlux();

staggeredLowerFaceScalarField phi("phi", fvMsh);

phi = ex::faceFlux(U);

// Adams-Bashforth source

staggeredScalarField H
(
    "H",
    fvMsh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE,
    true
);

H = ex::div(phi,U);
