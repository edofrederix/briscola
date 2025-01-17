#ifndef cellDataExchange_H
#define cellDataExchange_H

#include "dataExchange.H"
#include "meshField.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
class cellDataExchange
:
    public dataExchange<MeshType>
{
    // List of cell indices to obtain
    List<labelVector> indices_;

    // For each index, the receive number and index number
    List<Tuple2<label,label>> map_;

    // Neighboring processor numbers
    labelList neighbors_;

    // For each processor the cell indices to be sent
    List<List<labelVector>> sendCells_;

    // For each processor the cell indices to be received
    List<List<labelVector>> recvCells_;

private:

    // Disallow default bitwise assignment
    void operator=(const cellDataExchange&);

public:

    // Declare name of the class and its debug switch
    ClassName("cellDataExchange");

    // Construct from components
    cellDataExchange
    (
        const List<labelVector>&,
        const fvMesh&,
        const label l = 0,
        const label d = 0
    );

    // Copy constructor
    cellDataExchange(const cellDataExchange&);

    // Destructor
    virtual ~cellDataExchange();

    // (Re-)initialize from a local list of cell indices
    void init(const List<labelVector>&);

    // Access to the indices
    const List<labelVector>& indices() const { return indices_; }

    // Access to the map
    inline const List<Tuple2<label,label>>& map() const { return map_; }

    // Obtain the data

    template<class Type>
    List<Type> dataFunc(const meshField<Type,MeshType>& field) const;

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
