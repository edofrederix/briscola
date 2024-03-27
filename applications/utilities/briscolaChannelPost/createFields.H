staggeredScalarField stagField
(
    stag ? fieldName : "dummy",
    fvMsh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE,
    true
);

colocatedScalarField colField
(
    stag ? "dummy" : fieldName,
    fvMsh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE,
    true
);