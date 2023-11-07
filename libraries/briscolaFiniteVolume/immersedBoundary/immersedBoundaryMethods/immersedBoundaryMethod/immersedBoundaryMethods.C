#include "immersedBoundaryMethods.H"
#include "addToRunTimeSelectionTable.H"

#include "penalization.H"
#include "Fadlun.H"
#include "Vreman.H"
#include "Mittal.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeIBM(scalar,colocated);
makeIBM(scalar,staggered);
makeIBM(vector,colocated);
makeIBM(vector,staggered);

makeIBMType(penalization,scalar,colocated);
makeIBMType(penalization,scalar,staggered);
makeIBMType(penalization,vector,colocated);
makeIBMType(penalization,vector,staggered);

makeIBMType(Fadlun,scalar,colocated);
makeIBMType(Fadlun,scalar,staggered);
makeIBMType(Fadlun,vector,colocated);
makeIBMType(Fadlun,vector,staggered);

makeIBMType(Vreman,scalar,colocated);
makeIBMType(Vreman,scalar,staggered);
makeIBMType(Vreman,vector,colocated);
makeIBMType(Vreman,vector,staggered);

makeIBMType(Mittal,scalar,colocated);
makeIBMType(Mittal,scalar,staggered);
makeIBMType(Mittal,vector,colocated);
makeIBMType(Mittal,vector,staggered);

}

}

}
