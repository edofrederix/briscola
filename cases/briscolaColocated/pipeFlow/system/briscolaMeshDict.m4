FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      briscolaMeshDict;
}

vertices
(
    (-VARH -VARH 0)
    (VARH -VARH 0)
    (-VARH VARH 0)
    (VARH VARH 0)
    (-VARG -VARG 0)
    (VARG -VARG 0)
    (-VARG VARG 0)
    (VARG VARG 0)

    (-VARH -VARH VARL)
    (VARH -VARH VARL)
    (-VARH VARH VARL)
    (VARH VARH VARL)
    (-VARG -VARG VARL)
    (VARG -VARG VARL)
    (-VARG VARG VARL)
    (VARG VARG VARL)
);

bricks
(
    0
    {
        vertices    (0 8 2 10 1 9 3 11);
        N           (VARNX VARNY VARNZ);
        grading     geometric (1 1 1);
    }

    1
    {
        vertices    (4 12 0 8 5 13 1 9);
        N           (VARNX VARNR VARNZ);
        grading     geometric (1 VARGRAD 1);
    }

    2
    {
        vertices    (4 12 6 14 0 8 2 10);
        N           (VARNR VARNY VARNZ);
        grading     geometric (VARGRAD 1 1);
    }

    3
    {
        vertices    (1 9 3 11 5 13 7 15);
        N           (VARNR VARNY VARNZ);
        grading     geometric (VARGRADI 1 1);
    }

    4
    {
        vertices    (2 10 6 14 3 11 7 15);
        N           (VARNX VARNR VARNZ);
        grading     geometric (1 VARGRADI 1);
    }
);

edges
(
    (4 5)
    {
        type    arc;
        point   (0 -VARR 0);
    }

    (5 7)
    {
        type    arc;
        point   (VARR 0 0);
    }

    (7 6)
    {
        type    arc;
        point   (0 VARR 0);
    }

    (6 4)
    {
        type    arc;
        point   (-VARR 0 0);
    }

    (12 13)
    {
        type    arc;
        point   (0 -VARR VARL);
    }

    (13 15)
    {
        type    arc;
        point   (VARR 0 VARL);
    }

    (14 15)
    {
        type    arc;
        point   (0 VARR VARL);
    }

    (14 12)
    {
        type    arc;
        point   (-VARR 0 VARL);
    }
);

patches
(
    inlet
    {
        type        periodic;
        neighbor    outlet;

        faces
        (
            (0 1 2 3)
            (4 5 0 1)
            (4 0 6 2)
            (1 5 3 7)
            (2 3 6 7)
        );
    }

    outlet
    {
        type        periodic;
        neighbor    inlet;

        faces
        (
            (8 9 10 11)
            (12 13 8 9)
            (12 8 14 10)
            (9 13 15 11)
            (10 11 14 15)
        );
    }
);

decomposition
{
    type    manual;

    brickDecompositions
    (
        (1 1 4)
        (1 1 4)
        (1 1 4)
        (1 1 4)
        (1 1 4)
    );
}
