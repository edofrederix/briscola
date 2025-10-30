#include "RungeKuttaSchemes.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardEuler.H"
#include "forwardEuler.H"

#include "theta.H"
#include "CrankNicolson.H"

#include "RK3.H"
#include "RK4.H"
#include "ps3p5q.H"
#include "ps4p7q.H"

#include "DIRK2.H"
#include "DIRK3.H"
#include "DIRK4.H"
#include "Alexander3.H"
#include "Hairer5.H"
#include "midPoint.H"

#include "Ascher111.H"
#include "Ascher121.H"
#include "Ascher122.H"
#include "Ascher222.H"
#include "Ascher232.H"
#include "Ascher233.H"

#include "AB2.H"
#include "CNAB.H"
#include "CNABM.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

// First order schemes

makeRungeKuttaScheme(backwardEuler);
makeRungeKuttaScheme(forwardEuler);

// theta schemes

makeRungeKuttaScheme(theta);
makeRungeKuttaScheme(CrankNicolson);

// Explicit schemes

makeRungeKuttaScheme(RK3);
makeRungeKuttaScheme(RK4);
makeRungeKuttaScheme(ps3p5q);
makeRungeKuttaScheme(ps4p7q);

// Implicit schemes

makeRungeKuttaScheme(DIRK2);
makeRungeKuttaScheme(DIRK3);
makeRungeKuttaScheme(DIRK4);
makeRungeKuttaScheme(Alexander3);
makeRungeKuttaScheme(Hairer5);
makeRungeKuttaScheme(midPoint);

// IMEX schemes

makeRungeKuttaScheme(Ascher111);
makeRungeKuttaScheme(Ascher121);
makeRungeKuttaScheme(Ascher122);
makeRungeKuttaScheme(Ascher222);
makeRungeKuttaScheme(Ascher232);
makeRungeKuttaScheme(Ascher233);

// Multi-step schemes

makeRungeKuttaScheme(AB2);
makeRungeKuttaScheme(CNAB);
makeRungeKuttaScheme(CNABM);

}

}

}

}
