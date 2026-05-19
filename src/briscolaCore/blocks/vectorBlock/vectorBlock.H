#ifndef vectorBlock_H
#define vectorBlock_H

#include "block.H"
#include "vector.H"

#define TEMPLATE template<int P>
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "vectorBlock.C"
#endif

#endif
