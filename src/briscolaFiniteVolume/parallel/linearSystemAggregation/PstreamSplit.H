#ifndef PstreamSplit_H
#define PstreamSplit_H

#include "label.H"
#include "Tuple2.H"
#include "DynamicList.H"

#include "Pstream.H"
#include "PstreamGlobals.H"

namespace Foam
{

namespace PstreamGlobals
{

// List of split communicators
extern DynamicList<MPI_Comm> splitComms_;

// List to keep track of member procs
extern DynamicList<List<bool>> splitCommMembers_;

// Function to set a one-group split communicator and return its number
label split(const bool member);

// Return the split communicator from communicator number
MPI_Comm& getSplitComm(const label);

}

}

#endif
