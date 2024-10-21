#ifndef nearestPointInterpolator_H
#define nearestPointInterpolator_H

#include "pointInterpolator.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
class nearestPointInterpolator
:
    public pointInterpolator<MeshType>
{
    // Disallow default bitwise assignment
    void operator=(const nearestPointInterpolator&);

public:

    // Runtime type information
    TypeName("nearest");

    // Construct from components
    nearestPointInterpolator
    (
        const fvMesh&,
        const vectorList&,
        const bool global = true,
        const label l = 0,
        const label d = 0
    );

    // Copy constructors
    nearestPointInterpolator(const nearestPointInterpolator<MeshType>&);

    // Clone the interpolator
    virtual autoPtr<pointInterpolator<MeshType>> clone() const
    {
        return autoPtr<pointInterpolator<MeshType>>
        (
            new nearestPointInterpolator(*this)
        );
    }

    // Destructor
    virtual ~nearestPointInterpolator();

    // Interpolators

    #define INTERPFUNC(TYPE) \
    virtual List<TYPE> operator() (const meshField<TYPE,MeshType>&) const;

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
    List<Type> interp(const meshField<Type,MeshType>&) const;
};

}

}

}

#endif
