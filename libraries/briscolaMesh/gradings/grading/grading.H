#ifndef grading_H
#define grading_H

#include "label.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"
#include "labelTensor.H"

namespace Foam
{

namespace briscola
{

class brick;

class grading
{
protected:

    // Reference to the brick
    const brick& b_;

private:

    // Copy constructor
    grading(const grading&);

    // Disallow default bitwise assignment
    void operator=(const grading&);

public:

    // Runtime type information
    TypeName("grading");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        grading,
        dictionary,
        (
            const brick& b
        ),
        (b)
    );

    // Construct from components
    grading(const brick&);

    // Selector
    static autoPtr<grading> New(const brick&);

    // Destructor
    virtual ~grading();

    // Mappig function from uniform to graded for a list of normalized
    // coordinates and for the given edge
    virtual tmp<scalarField> operator()
    (
        const scalarField&,
        const label
    ) const = 0;

    // Transform the grading given a transformation tensor
    virtual void transform(const labelTensor) = 0;
};

}

}

#endif
