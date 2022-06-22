#ifndef meshObject_H
#define meshObject_H

#include "label.H"

namespace Foam
{

namespace briscola
{

class edge;
class face;
class brick;
class geometry;

template<class ParentType>
class meshObject
{
protected:

    // Reference to the parent object
    const ParentType& parent_;

    // Number of this object in its parent
    const label num_;

public:

    // Construct from components
    meshObject(const ParentType&, const label);

    // Destructor
    virtual ~meshObject();

    // Access to the parent
    inline const ParentType& parent() const
    {
        return parent_;
    }

    // Return the number of this object in its parent
    inline label num() const
    {
        return num_;
    }
};

}

}

#endif
