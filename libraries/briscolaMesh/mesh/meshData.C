#include "meshData.H"

namespace Foam
{

namespace briscola
{

void meshData::checkData() const
{
    forAll(vertexData_, i)
    {
        if (vertexData_[i].size() != 3)
        {
            FatalErrorInFunction
                << "Incorrect number of components specified for vertex "
                << i << ": " << vertexData_[i] << endl
                << exit(FatalError);
        }
    }

    forAll(brickData_.sortedToc(), i)
    {
        const label key = brickData_.sortedToc()[i];

        if (key != i)
        {
            FatalErrorInFunction
                << "Incorrect labeling of brick dictionaries. The "
                << i << "th brick has label " << key << endl
                << exit(FatalError);
        }

        const labelBlock vertices(brickData_[i].lookup("vertices"));

        if (vertices.size() != 8)
        {
            FatalErrorInFunction
                << "Incorrect number of vertices specified for brick "
                << i << ": " << vertices << endl
                << exit(FatalError);
        }

        forAllBlockLinear(vertices, l)
        {
            if (vertices(l) < 0 || vertices(l) >= vertexData_.size())
            {
                FatalErrorInFunction
                    << "Incorrect vertex label found in brick "
                    << i << ": " << vertices << endl
                    << exit(FatalError);
            }
        }
    }

    forAllConstIter(patchTable, patchData_, iter)
    {
        const dictionary& patchDict = iter();

        const word name(iter.key());

        const List<labelList> vertexNums(patchDict.lookup("faces"));

        forAll(vertexNums, i)
        {
            forAll(vertexNums[i], j)
            {
                if (vertexNums[i][j] < 0 || vertexNums[i][j] >= vertexData_.size())
                {
                    FatalErrorInFunction
                        << "Incorrect vertex label found in face "
                        << i << " of patch " << name << endl
                        << exit(FatalError);
                }
            }
        }
    }
}

meshData::meshData(const IOdictionary& dict)
:
    dict_(dict),
    vertexData_(dict.lookup("vertices")),
    brickData_(dict.lookup("bricks")),
    patchData_(dict.lookup("patches")),
    edgeData_(dict.lookup("edges"))
{
    checkData();
}

meshData::~meshData()
{}

}

}
