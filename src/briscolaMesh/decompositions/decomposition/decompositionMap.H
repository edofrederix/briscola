#ifndef decompositionMap_H
#define decompositionMap_H

#include "labelBlock.H"
#include "List.H"
#include "labelVector.H"

namespace Foam
{

namespace briscola
{

class decompositionMap
:
    public labelBlock
{
    // Brick legend
    List<labelVector> legend_;

    // Reset and calculate legend from map
    void calcLegend();

public:

    // Null constructor
    decompositionMap();

    // Construct from map
    decompositionMap(const labelBlock&);

    // Destructor
    virtual ~decompositionMap();

    // Access to the brick legend
    const List<labelVector>& legend() const { return legend_; }

    // Set data from map and calculate legend
    void setData(const labelBlock&);
};

}

}

#endif
