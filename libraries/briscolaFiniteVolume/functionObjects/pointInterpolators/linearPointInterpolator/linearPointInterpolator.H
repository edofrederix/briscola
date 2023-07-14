#ifndef linearPointInterpolator_H
#define linearPointInterpolator_H

#include "pointInterpolator.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class linearPointInterpolator
:
    public pointInterpolator
{
protected:

    // Interpolation weights
    List<vertexScalar> weights_;

private:

    // Disallow default bitwise assignment
    void operator=(const linearPointInterpolator&);

public:

    // Runtime type information
    TypeName("linear");

    // Construct from components
    linearPointInterpolator(const fvMesh&, const vectorList&);

    // Copy constructors
    linearPointInterpolator(const linearPointInterpolator&);

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
