#ifndef pointDataExchange_H
#define pointDataExchange_H

#include "dataExchange.H"
#include "meshField.H"
#include "pointInterpolator.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
class pointDataExchange
:
    public dataExchange<MeshType>
{
    // Point interpolator
    autoPtr<pointInterpolator<MeshType>> interpPtr_;

private:

    // Disallow default bitwise assignment
    void operator=(const pointDataExchange&);

public:

    // Declare name of the class and its debug switch
    ClassName("pointDataExchange");

    // Construct from components
    pointDataExchange
    (
        const List<vector>&,
        const fvMesh&,
        const label l = 0,
        const label d = 0
    );

    // Copy constructor
    pointDataExchange(const pointDataExchange&);

    // Destructor
    virtual ~pointDataExchange();

    // (Re-)initialize from a local list of points
    void init(const List<vector>&);

    // Access to the interpolator
    const pointInterpolator<MeshType>& interp() const { return interpPtr_; }

    // Obtain the data

    template<class Type>
    List<Type> dataFunc(const meshField<Type,MeshType>& field);

    #define DATAFUNC(TYPE)                                                     \
    virtual List<TYPE> operator() (const meshField<TYPE,MeshType>& field)      \
    {                                                                          \
        return this->dataFunc(field);                                          \
    }

    DATAFUNC(scalar);
    DATAFUNC(vector);
    DATAFUNC(tensor);
    DATAFUNC(sphericalTensor);
    DATAFUNC(symmTensor);
    DATAFUNC(diagTensor);
    DATAFUNC(faceScalar);
    DATAFUNC(faceVector);

    #undef DATAFUNC
};

}

}

}

#endif
