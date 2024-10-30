FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      U;
}

boundaryConditions
{
    walls
    {
        type        noSlip;
    }

    inlet
    {
        type        periodic;
    }

    outlet
    {
        type        periodic;
    }

    ib
    {
        type        VARIBMBC;
        values      ((0 0 0));
    }
}
