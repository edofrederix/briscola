#include "decompositionMap.H"
#include "brick.H"
#include "Tuple2.H"

namespace Foam
{

namespace briscola
{

decompositionMap::decompositionMap()
:
    labelBlock(),
    legend_()
{}

decompositionMap::decompositionMap(const labelBlock& map)
:
    labelBlock(map),
    legend_()
{
    calcLegend();
}

decompositionMap::~decompositionMap()
{}

void decompositionMap::calcLegend()
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

void decompositionMap::setData(const labelBlock& map)
{
    *this = map;
    calcLegend();
}

}

}
