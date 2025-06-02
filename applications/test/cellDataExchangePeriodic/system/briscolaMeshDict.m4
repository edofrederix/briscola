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
        N           (12 12 12);
    }
);

edges
();

patches
(
    aft
    {
        type        periodic;
        neighbor    fore;

        faces
        (
            (0 1 2 3)
        );
    }

    fore
    {
        type        periodic;
        neighbor    aft;

        faces
        (
            (4 5 6 7)
        );
    }

    bottom
    {
        type        periodic;
        neighbor    top;

        faces
        (
            (0 1 4 5)
        );
    }

    top
    {
        type        periodic;
        neighbor    bottom;

        faces
        (
            (2 3 6 7)
        );
    }

    left
    {
        type        periodic;
        neighbor    right;

        faces
        (
            (0 4 2 6)
        );
    }

    right
    {
        type        periodic;
        neighbor    left;

        faces
        (
            (1 5 3 7)
        );
    }
);

decomposition
{
    type    manual;

    brickDecompositions
    (
        (VARNPROC VARNPROC VARNPROC)
    );
}
