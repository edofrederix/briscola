FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      briscolaMeshDict;
}

vertices
(
    (0 0 0)
    (1 0 0)
    (0 1 0)
    (1 1 0)

    (0 0 1)
    (1 0 1)
    (0 1 1)
    (1 1 1)
);

bricks
(
    0
    {
        vertices    (((0 4)(2 6))((1 5)(3 7)));
        N           (VARN VARN VARN);
    }
);

edges
();

patches
(
    walls
    {
        type    patch;

        faces
        (
            (0 1 4 5)
            (2 3 6 7)
        );
    }

    inlet
    {
        type    patch;

        faces
        (
            (0 2 4 6)
        );
    }

    outlet
    {
        type    patch;

        faces
        (
            (1 3 5 7)
        );
    }

    frontAndBack
    {
        type    patch;

        faces
        (
            (0 1 2 3)
            (4 5 6 7)
        );
    }

);

decomposition
{
    type    manual;

    brickDecompositions
    (
        (VARNP VARNP VARNP)
    );
}
