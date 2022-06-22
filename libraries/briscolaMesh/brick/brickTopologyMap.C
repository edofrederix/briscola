#include "brickTopologyMap.H"
#include "brick.H"
#include "Tuple2.H"

namespace Foam
{

namespace briscola
{

brickTopologyMap::brickTopologyMap()
:
    map_(),
    legend_()
{}

brickTopologyMap::brickTopologyMap(const labelBlock& map)
:
    map_(map)
{
    calcLegend();
}

brickTopologyMap::~brickTopologyMap()
{}

void brickTopologyMap::calcLegend()
{
    legend_.clear();

    forAllBlock(map_, i, j, k)
    if (map_(i,j,k) != -1)
    {
        label v = map_(i,j,k);

        if (v >= legend_.size())
        {
            legend_.setSize(v+1);
        }

        legend_[v] = labelVector(i,j,k);
    }
}

void brickTopologyMap::setData(const labelBlock& map)
{
    map_ = map;
    calcLegend();
}

}

}
