FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      briscolaMeshDict;
}

vertices
(
    (0 -1.2 0)
    (8 -1.2 0)
    (0  1.2 0)
    (8  1.2 0)

    (0 -1.2 0.1)
    (8 -1.2 0.1)
    (0  1.2 0.1)
    (8  1.2 0.1)
);

bricks
(
    0
    {
        vertices    (((0 4)(2 6))((1 5)(3 7)));
        N           (VARMESHX VARMESHY 1);
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
        type    periodic;
        neighbor outlet;

        faces
        (
            (0 2 4 6)
        );
    }

    outlet
    {
        type    periodic;
        neighbor inlet;

        faces
        (
            (1 3 5 7)
        );
    }

    empties
    {
        type    empty;

        faces
        (
            (0 1 2 3)
            (4 5 6 7)
        );
    }
);

immersedBoundaries
{
    ib
    {
        cylinder
        {
            type        cylinder;
            start       (0 0 0);
            end         (8 0 0);
            radius      1;
            inverted    true;
        }
    }
}

decomposition
{
    type    manual;

    brickDecompositions
    (
        (VARNPROCX VARNPROCY 1)
    );
}
