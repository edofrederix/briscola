FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      briscolaMeshDict;
}

vertices
(
    (-VARRSQRT -VARRSQRT -VARW)
    ( VARRSQRT -VARRSQRT -VARW)
    (-VARRSQRT  VARRSQRT -VARW)
    ( VARRSQRT  VARRSQRT -VARW)

    (-VARR2SQRT -VARR2SQRT -VARW)
    ( VARR2SQRT -VARR2SQRT -VARW)
    (-VARR2SQRT  VARR2SQRT -VARW)
    ( VARR2SQRT  VARR2SQRT -VARW)

    (-VARLI     -VARH -VARW)
    (-VARR2SQRT -VARH -VARW)
    ( VARR2SQRT -VARH -VARW)
    ( VARLO     -VARH -VARW)

    (-VARLI     -VARR2SQRT -VARW)
    ( VARLO     -VARR2SQRT -VARW)

    (-VARLI      VARR2SQRT -VARW)
    ( VARLO      VARR2SQRT -VARW)

    (-VARLI      VARH -VARW)
    (-VARR2SQRT  VARH -VARW)
    ( VARR2SQRT  VARH -VARW)
    ( VARLO      VARH -VARW)


    (-VARRSQRT -VARRSQRT VARW)
    ( VARRSQRT -VARRSQRT VARW)
    (-VARRSQRT  VARRSQRT VARW)
    ( VARRSQRT  VARRSQRT VARW)

    (-VARR2SQRT -VARR2SQRT VARW)
    ( VARR2SQRT -VARR2SQRT VARW)
    (-VARR2SQRT  VARR2SQRT VARW)
    ( VARR2SQRT  VARR2SQRT VARW)

    (-VARLI     -VARH VARW)
    (-VARR2SQRT -VARH VARW)
    ( VARR2SQRT -VARH VARW)
    ( VARLO     -VARH VARW)

    (-VARLI     -VARR2SQRT VARW)
    ( VARLO     -VARR2SQRT VARW)

    (-VARLI      VARR2SQRT VARW)
    ( VARLO      VARR2SQRT VARW)

    (-VARLI      VARH VARW)
    (-VARR2SQRT  VARH VARW)
    ( VARR2SQRT  VARH VARW)
    ( VARLO      VARH VARW)
);

bricks
(
    0
    {
        vertices    (4 24 6 26 0 20 2 22);
        N           (VARNR VARNQ 1);
    }

    1
    {
        vertices    (1 21 3 23 5 25 7 27);
        N           (VARNR VARNQ 1);
    }

    2
    {
        vertices    (4 24 0 20 5 25 1 21);
        N           (VARNQ VARNR 1);
    }

    3
    {
        vertices    (2 22 6 26 3 23 7 27);
        N           (VARNQ VARNR 1);
    }

    4
    {
        vertices    (8 28 12 32 9 29 4 24);
        N           (VARNI VARNS 1);
        grading     geometric (VARGIi VARGSi 1);
    }

    5
    {
        vertices    (9 29 4 24 10 30 5 25);
        N           (VARNQ VARNS 1);
        grading     geometric (1 VARGSi 1);
    }

    6
    {
        vertices    (10 30 5 25 11 31 13 33);
        N           (VARNO VARNS 1);
        grading     geometric (VARGO VARGSi 1);
    }

    7
    {
        vertices    (12 32 14 34 4 24 6 26);
        N           (VARNI VARNQ 1);
        grading     geometric (VARGIi 1 1);
    }

    8
    {
        vertices    (5 25 7 27 13 33 15 35);
        N           (VARNO VARNQ 1);
        grading     geometric (VARGO 1 1);
    }

    9
    {
        vertices    (14 34 16 36 6 26 17 37);
        N           (VARNI VARNS 1);
        grading     geometric (VARGIi VARGS 1);
    }

    10
    {
        vertices    (6 26 17 37 7 27 18 38);
        N           (VARNQ VARNS 1);
        grading     geometric (1 VARGS 1);
    }

    11
    {
        vertices    (7 27 18 38 15 35 19 39);
        N           (VARNO VARNS 1);
        grading     geometric (VARGO VARGS 1);
    }
);

edges
(
    (0 1)
    {
        type    arc;
        point   (0 -VARR -VARW);
    }

    (1 3)
    {
        type    arc;
        point   (VARR 0 -VARW);
    }

    (3 2)
    {
        type    arc;
        point   (0 VARR -VARW);
    }

    (2 0)
    {
        type    arc;
        point   (-VARR 0 -VARW);
    }


    (20 21)
    {
        type    arc;
        point   (0 -VARR VARW);
    }

    (21 23)
    {
        type    arc;
        point   (VARR 0 VARW);
    }

    (23 22)
    {
        type    arc;
        point   (0 VARR VARW);
    }

    (22 20)
    {
        type    arc;
        point   (-VARR 0 VARW);
    }


    (4 5)
    {
        type    arc;
        point   (0 -VARR2 -VARW);
    }

    (5 7)
    {
        type    arc;
        point   (VARR2 0 -VARW);
    }

    (7 6)
    {
        type    arc;
        point   (0 VARR2 -VARW);
    }

    (6 4)
    {
        type    arc;
        point   (-VARR2 0 -VARW);
    }


    (24 25)
    {
        type    arc;
        point   (0 -VARR2 VARW);
    }

    (25 27)
    {
        type    arc;
        point   (VARR2 0 VARW);
    }

    (27 26)
    {
        type    arc;
        point   (0 VARR2 VARW);
    }

    (26 24)
    {
        type    arc;
        point   (-VARR2 0 VARW);
    }
);

patches
(
    inlet
    {
        type        patch;

        faces
        (
            (8 28 12 32)
            (12 32 14 34)
            (14 34 16 36)
        );
    }

    outlet
    {
        type        patch;

        faces
        (
            (11 31 13 33)
            (13 33 15 35)
            (15 35 19 39)
        );
    }

    wall
    {
        type        patch;

        faces
        (
            (0 20 1 21)
            (1 21 3 23)
            (2 22 3 23)
            (0 20 2 22)
        );
    }

    bottom
    {
        type        patch;

        faces
        (
            (8 28 9 29)
            (9 29 10 30)
            (10 30 11 31)
        );
    }

    top
    {
        type        patch;

        faces
        (
            (16 36 17 37)
            (17 37 18 38)
            (18 38 19 39)
        );
    }

    empties
    {
        type        empty;

        faces
        (
            (8 12 9 4)
            (9 4 10 5)
            (10 5 11 13)
            (12 14 4 6)
            (4 6 0 2)
            (4 0 5 1)
            (1 3 5 7)
            (5 7 13 15)
            (2 6 3 7)
            (14 16 6 17)
            (6 17 7 18)
            (7 18 15 19)
            (28 32 29 24)
            (29 24 30 25)
            (30 25 31 33)
            (32 34 24 26)
            (24 26 20 22)
            (24 20 25 21)
            (21 23 25 27)
            (22 26 23 27)
            (25 27 33 35)
            (34 36 26 37)
            (26 37 27 38)
            (27 38 35 39)
        );
    }
);

decomposition
{
    type    manual;

    brickDecompositions
    (
        (1 1 1)
        (1 1 1)
        (1 1 1)
        (1 1 1)
        (1 1 1)
        (1 1 1)
        (2 1 1)
        (1 1 1)
        (2 1 1)
        (1 1 1)
        (1 1 1)
        (2 1 1)
    );
}
