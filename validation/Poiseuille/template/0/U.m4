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
        setAverage  true;
        average     (1 0 0);
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
        value       (0 0 0);
    }
}
