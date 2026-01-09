#ifndef decompositionSlice_H
#define decompositionSlice_H

#include "label.H"
#include "labelVector.H"

namespace Foam
{

namespace briscola
{

class decomposition;

class decompositionSlice
{
    // Processor number in the first brick
    const label procNum0_;

    // Processor number in the second brick
    const label procNum1_;

    // Index to the start point of the slice in the first brick
    const labelVector start0_;

    // Index to the start point of the slice in the second brick
    const labelVector start1_;

    // Size of the slice in the first brick
    const labelVector N0_;

    // Size of the slice in the second brick
    const labelVector N1_;

    // Disallow default bitwise assignment
    void operator=(const decompositionSlice&);

public:

    // Construct from components
    decompositionSlice
    (
        const label,
        const label,
        const labelVector&,
        const labelVector&,
        const labelVector&,
        const labelVector&
    );

    // Copy constructor
    decompositionSlice(const decompositionSlice&);

    // Destructor
    virtual ~decompositionSlice();

    // Return the processor numbers
    inline label procNum0() const { return procNum0_; }
    inline label procNum1() const { return procNum1_; }

    // Return the start point of the slice
    inline labelVector start0() const { return start0_; }
    inline labelVector start1() const { return start1_; }

    // Return the size of the slice
    inline labelVector N0() const { return N0_; }
    inline labelVector N1() const { return N1_; }
};

}

}

#endif
