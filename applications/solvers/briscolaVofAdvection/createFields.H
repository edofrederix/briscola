// Velocity

colocatedVectorField U
(
    "U",
    fvMsh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE,
    true,
    true
);

U = Zero;

// Flux

colocatedFaceScalarField phi
(
    "phi",
    fvMsh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE,
    true
);

phi = Zero;

//  Normal

colocatedVectorField n
(
    "n",
    fvMsh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE,
    true
);
