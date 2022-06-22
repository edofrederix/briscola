#ifndef meshType_H
#define meshType_H

#include "labelVector.H"
#include "vector.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class meshType
{
public:

    // Number of directions
    static const int numberOfDirections;

    // Direction names
    static const char* const directionNames[];

    // Padding is an array of labelVectors, giving for each mesh direction the
    // number of cells that must be added in each spatial dimension. Padding
    // does not include ghost cell layers.
    static const labelVector padding[];

    // Shift is an array of vectors, giving for each mesh direction the shift of
    // cells w.r.t. the base mesh
    static const vector shift[];

    // Mesh type name
    static const char* const typeName;
};

}

}

}

#endif
