#ifndef nearestPointInterpolator_H
#define nearestPointInterpolator_H

#include "pointInterpolator.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class nearestPointInterpolator
:
    public pointInterpolator
{
    // Disallow default bitwise assignment
    void operator=(const nearestPointInterpolator&);

public:

    // Runtime type information
    TypeName("nearest");

    // Construct from components
    nearestPointInterpolator(const fvMesh&, const vectorList&);

    // Copy constructors
    nearestPointInterpolator(const nearestPointInterpolator&);

    // Interpolators

    #define INTERPFUNC(TYPE) \
    virtual List<TYPE> operator() (const meshField<TYPE,colocated>&);

    INTERPFUNC(scalar);
    INTERPFUNC(vector);
    INTERPFUNC(tensor);
    INTERPFUNC(sphericalTensor);
    INTERPFUNC(symmTensor);
    INTERPFUNC(diagTensor);
    INTERPFUNC(faceScalar);
    INTERPFUNC(faceVector);

    #undef INTERPFUNC

    // Colocated interpolator
    template<class Type>
    List<Type> interp(const meshField<Type,colocated>&);
};

}

}

}

#endif
