#ifndef dataExchange_H
#define dataExchange_H

#include "fvMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
class dataExchange
{
protected:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Mesh level
    const label l_;

    // Mesh direction
    const label d_;

private:

    // Disallow default bitwise assignment
    void operator=(const dataExchange&);

public:

    // Declare name of the class and its debug switch
    ClassName("dataExchange");

    // Construct from mesh and mesh direction
    dataExchange(const fvMesh&, const label l = 0, const label d = 0);

    // Copy constructor
    dataExchange(const dataExchange&);

    // Destructor
    virtual ~dataExchange();

    // Return access to the mesh
    inline const fvMesh& fvMsh() const
    {
        return fvMsh_;
    }

    // Return the level number
    inline label l() const { return l_; }

    // Return the direction number
    inline label d() const { return d_; }

    // Obtain the data

    #define DATAFUNC(TYPE) \
    virtual List<TYPE> operator() (const meshField<TYPE,MeshType>&) = 0;

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
