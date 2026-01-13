#ifndef CNABM_H
#define CNABM_H

#include "RungeKuttaScheme.H"

/**

@brief Modified Crank-Nicolson/Adams-Bashforth (CNAB) IMEX scheme

Crank-Nicolson/Adams-Bashforth (CNAB) IMEX scheme, see Eq. (2) of Ascher, U. M.,
Ruuth, S. J., & Wetton, B. T. (1995). Implicit-explicit methods for
time-dependent partial differential equations. SIAM Journal on Numerical
Analysis, 32(3), 797-823.

*/

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class CNABM
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("CNABM");

    // Construct from mesh
    CNABM(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "CNABM.C"
#endif

#endif
