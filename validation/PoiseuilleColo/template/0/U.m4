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

    inlet_mapped
    {
        type        mapped;
        setAverages (true);
        averages    ((1 0 0));
    }

    inlet
    {
        type        periodic;
    }

    outlet_mapped
    {
        type        outflow;
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
