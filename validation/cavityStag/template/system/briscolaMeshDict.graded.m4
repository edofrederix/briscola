FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      briscolaMeshDict;
}

vertices
(
    (0   0   0)         // 0
    (0.5 0   0)         // 1
    (1   0   0)         // 2

    (0   0.5 0)         // 3
    (0.5 0.5 0)         // 4
    (1   0.5 0)         // 5

    (0   1   0)         // 6
    (0.5 1   0)         // 7
    (1   1   0)         // 8

    (0   0   0.1)       // 9
    (0.5 0   0.1)       // 10
    (1   0   0.1)       // 11

    (0   0.5 0.1)       // 12
    (0.5 0.5 0.1)       // 13
    (1   0.5 0.1)       // 14

    (0   1   0.1)       // 15
    (0.5 1   0.1)       // 16
    (1   1   0.1)       // 17
);

bricks
(
    0
    {
        vertices    (((0 9)(3 12))((1 10)(4 13)));
        N           (48 48 1);
        grading     geometric (4 4 1);
    }

    1
    {
        vertices    (((3 12)(6 15))((4 13)(7 16)));
        N           (48 48 1);
        grading     geometric (4 0.25 1);
    }

    2
    {
        vertices    (((1 10)(4 13))((2 11)(5 14)));
        N           (48 48 1);
        grading     geometric (0.25 4 1);
    }

    3
    {
        vertices    (((4 13)(7 16))((5 14)(8 17)));
        N           (48 48 1);
        grading     geometric (0.25 0.25 1);
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
            (0 9 3 12)
            (3 12 6 15)
            (2 11 5 14)
            (5 14 8 17)
            (0 9 1 10)
            (1 10 2 11)
        );
    }

    movingWall
    {
        type    patch;

        faces
        (
            (6 15 7 16)
            (7 16 8 17)
        );
    }
    
    empties
    {
        type    empty;

        faces
        (
            (9 10 12 13)
            (12 13 15 16)
            (10 11 13 14)
            (13 14 16 17)
            (0 1 3 4)
            (1 2 4 5)
            (3 4 6 7)
            (4 5 7 8)
        );
    }

);

decomposition
{
    type    manual;

    brickDecompositions
    (
        (VARNPROCX VARNPROCY 1)
        (VARNPROCX VARNPROCY 1)
        (VARNPROCX VARNPROCY 1)
        (VARNPROCX VARNPROCY 1)
    );
}
