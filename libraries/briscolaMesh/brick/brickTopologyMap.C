#include "brickTopologyMap.H"
#include "brick.H"
#include "Tuple2.H"

namespace Foam
{

namespace briscola
{

brickTopologyMap::brickTopologyMap()
:
    labelBlock(),
    legend_()
{}

brickTopologyMap::brickTopologyMap(const labelBlock& map)
:
    labelBlock(map)
{
    calcLegend();
}

brickTopologyMap::~brickTopologyMap()
{}

void brickTopologyMap::calcLegend()
{
    legend_.clear();

    labelBlock& map = *this;

    forAllBlock(map, i, j, k)
    if (map(i,j,k) != -1)
    {
        label v = map(i,j,k);

        if (v >= legend_.size())
        {
            legend_.setSize(v+1);
        }

        legend_[v] = labelVector(i,j,k);
    }
}

void brickTopologyMap::setData(const labelBlock& map)
{
    *this = map;
    calcLegend();
}

}

}
