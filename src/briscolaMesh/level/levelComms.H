#ifndef levelComms_H
#define levelComms_H

#include "labelList.H"
#include "labelBlock.H"

// Forward declaration of MPI communicator
struct ompi_communicator_t;
typedef ompi_communicator_t* MPI_Comm;

namespace Foam
{

namespace briscola
{

class level;

class levelComms
{
    // Level to which these boundaries belong
    const level& lvl_;

    // Level communication member data

        // Level communicator number
        label comm_;

        // Neighbors
        labelList recvNeighbors_;
        labelList sendNeighbors_;

        // Neighbor maps
        List<labelList> recvMap_;
        List<labelList> sendMap_;

    // Agglomerate communicator number which is defined on the agglomerated
    // level
    label agg_;

    // Access to the level communicator
    MPI_Comm& communicator();

    // Initialize the level communicator
    void initLevelComm();

    // Initialize the agglomerate communicator
    void initAgglomerateComm();

public:

    // Construct from level
    levelComms(const level&);

    // Copy constructor
    levelComms(const levelComms&);

    // Copy constructor with alternative level reference
    levelComms(const levelComms&, const level&);

    // Destructor
    ~levelComms();

    // Level
    inline const level& lvl() const { return lvl_; }

    // Level communicator number
    label comm() const { return comm_; }

    // Implicit conversion to the communicator number
    operator label() const { return comm_; }

    // Return const reference to the level communicator
    const MPI_Comm& communicator() const;

    // Return const reference to the agglomerate communicator
    const MPI_Comm& aggCommunicator() const;

    // Receive neighbors
    const labelList& recvNeighbors() const { return recvNeighbors_; }

    // Send neighbors
    const labelList& sendNeighbors() const { return sendNeighbors_; }

    // Receive neighbors
    const List<labelList>& recvMap() const { return recvMap_; }

    // Send neighbors
    const List<labelList>& sendMap() const { return sendMap_; }

    // Agglomerate communicator number
    label agg() const { return agg_; }
};

}

}

#endif
