#ifndef PstreamGlobalsLsa_H
#define PstreamGlobalsLsa_H

#include "label.H"
#include "Tuple2.H"
#include "DynamicList.H"

#include "Pstream.H"
#include "PstreamGlobals.H"

namespace Foam
{

namespace PstreamGlobals
{

// List to keep track of created communicators and their indices
extern DynamicList<Tuple2<label,label>> lsaCommsData_;
extern DynamicList<Tuple2<label,label>> lsaMasterCommsData_;

// Function to set an aggregation communicator and return its number. Optionally
// set a communicator between masters.
label lsaSetComms
(
    const label nParts,
    const label color,
    const label key,
    const bool masters = false
);

// Return the aggregation communicator from communicator number
MPI_Comm& lsaGetComm(const label commNum_);

}

}

#endif
