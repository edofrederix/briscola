#ifndef linearPointInterpolator_H
#define linearPointInterpolator_H

#include "pointInterpolator.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
class linearPointInterpolator
:
    public pointInterpolator<MeshType>
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
    linearPointInterpolator
    (
        const fvMesh&,
        const vectorList&,
        const bool global = true,
        const label l = 0,
        const label d = 0
    );

    // Copy constructors
    linearPointInterpolator(const linearPointInterpolator<MeshType>&);

    // Clone the interpolator
    virtual autoPtr<pointInterpolator<MeshType>> clone() const
    {
        return autoPtr<pointInterpolator<MeshType>>
        (
            new linearPointInterpolator(*this)
        );
    }

    // Destructor
    virtual ~linearPointInterpolator();

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
