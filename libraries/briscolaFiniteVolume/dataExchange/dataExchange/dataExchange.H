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

    // Mesh direction
    const label d_;

    // List of processors to send to
    labelList sendProcNums_;

    // List of processors ro receive from
    labelList recvProcNums_;

    // List of send counts per processor
    labelList sendCounts_;

    // List of receive counts per processor
    labelList recvCounts_;

    // List of patch pointers corresponding to send processors
    List<const partPatch*> sendPatchPtrs_;

    // List of patch pointers corresponding to receive processors
    List<const partPatch*> recvPatchPtrs_;

    // Given a local list of processor numbers to talk to and receive counts,
    // setup send/recv lists
    void setupComms(const labelList&, const labelList&);

private:

    // Disallow default bitwise assignment
    void operator=(const dataExchange&);

public:

    // Declare name of the class and its debug switch
    ClassName("dataExchange");

    // Construct from mesh and mesh direction
    dataExchange(const fvMesh&, const label d = 0);

    // Copy constructor
    dataExchange(const dataExchange&);

    // Destructor
    virtual ~dataExchange();

    // Return access to the mesh
    inline const fvMesh& fvMsh() const
    {
        return fvMsh_;
    }

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
