#ifndef brickTopologyMap_H
#define brickTopologyMap_H

#include "labelBlock.H"
#include "List.H"
#include "labelVector.H"

namespace Foam
{

namespace briscola
{

class brickTopologyMap
{
    // Brick map
    labelBlock map_;

    // Brick legend
    List<labelVector> legend_;

    // Reset and calculate legend from map
    void calcLegend();

public:

    // Null constructor
    brickTopologyMap();

    // Construct from map
    brickTopologyMap(const labelBlock&);

    // Destructor
    virtual ~brickTopologyMap();

    // Access to the brick map
    const labelBlock& map() const { return map_; }

    // Access to the brick legend
    const List<labelVector>& legend() const { return legend_; }

    // Set data from map and calculate legend
    void setData(const labelBlock&);
};

}

}

#endif
