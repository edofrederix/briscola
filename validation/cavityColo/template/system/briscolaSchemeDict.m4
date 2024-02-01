FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      briscolaSchemesDict;
}

ddtSchemes
{
    ddt(U)
    {
        type    Euler;
    }
}

laplacianSchemes
{
    laplacian(p)
    {
        type    VARLSCHEME;
    }

    laplacian(const,U)
    {
        type    VARLSCHEME;
    }
}

divergenceSchemes
{
    div(phi,U)
    {
        type    VARDSCHEME;
    }
}

gradientSchemes
{
    grad(p)
    {
        type    midPointGauss;
    }
}

faceGradientSchemes
{
    faceGrad(p)
    {
        type    linear;
    }
}

faceFluxSchemes
{
    flux(U)
    {
        type    midPoint;
    }
}
