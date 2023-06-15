#ifndef shape_H
#define shape_H

#include "vector.H"

namespace Foam
{

namespace briscola
{

namespace ibm
{

class shape
{

public:

    // Empty constructor
    shape();

    // Destructor
    virtual ~shape();

    // Check if point is inside shape
    virtual bool isInside(vector) = 0;

};

} // end namespace ibm

} // end namespace briscola

} // end namespace Foam

#endif
