#ifndef meshData_H
#define meshData_H

#include "IOdictionary.H"
#include "labelList.H"
#include "vectorList.H"
#include "block.H"
#include "edgeKey.H"

namespace Foam
{

namespace briscola
{

typedef HashTable<dictionary, word, word::hash> patchTable;
typedef HashTable<dictionary, edgeKey, edgeKey::hash> edgeTable;
typedef HashTable<dictionary, label, Hash<label>> brickTable;

class meshData
{
    // Reference to the mesh dictionary
    const IOdictionary& dict_;

    // Vertex data
    vectorList vertexData_;

    // Brick data
    brickTable brickData_;

    // Patch data
    patchTable patchData_;

    // Edge data
    edgeTable edgeData_;

    // Check the raw input data
    void checkData() const;

public:

    // Construct from IOdictionary
    meshData(const IOdictionary&);

    // Copy constructor
    meshData(const meshData& data);

    // Destructor
    ~meshData();

    // Acces to the mesh dictionary
    inline const IOdictionary& dict() const
    {
        return dict_;
    }

    // Access to the vertex data
    inline const vectorList& vertexData() const
    {
        return vertexData_;
    }

    // Access to the brick data
    inline const brickTable& brickData() const
    {
        return brickData_;
    }

    // Access to the patch data
    inline const patchTable& patchData() const
    {
        return patchData_;
    }

    // Access to the edge data
    inline const edgeTable& edgeData() const
    {
        return edgeData_;
    }
};

}

}

#endif
