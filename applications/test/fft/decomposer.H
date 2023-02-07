#ifndef DECOMPOSER_H
#define DECOMPOSER_H

#include "vector.H"
#include "block.H"

namespace Foam
{

namespace briscola
{

class Decomposer
{

	//Global Size

	labelVector N_;

    // Initial decomposition

    labelVector I_;

    public:

    // X- Y- and Z-pencil decompositions

    labelVector X_;
    labelVector Y_;
    labelVector Z_;

    // Number of cells in each decomposition

    List<labelVector> Ni_;
    List<labelVector> Nx_;
    List<labelVector> Ny_;
    List<labelVector> Nz_;

    // Start cell indices of each processor in each decomposition

    List<labelVector> si_;
    List<labelVector> sx_;
    List<labelVector> sy_;
    List<labelVector> sz_;

    public:

	Decomposer(labelVector n, labelVector i)
    {
	    N_ = n;
        I_ = i;

        Ni_ = List<labelVector>(Pstream::nProcs(), Zero);
        Nx_ = List<labelVector>(Pstream::nProcs(), Zero);
        Ny_ = List<labelVector>(Pstream::nProcs(), Zero);
        Nz_ = List<labelVector>(Pstream::nProcs(), Zero);

        si_ = List<labelVector>(Pstream::nProcs(), Zero);
        sx_ = List<labelVector>(Pstream::nProcs(), Zero);
        sy_ = List<labelVector>(Pstream::nProcs(), Zero);
        sz_ = List<labelVector>(Pstream::nProcs(), Zero);

        decompInit();
	}

    // Decomposition functions

    void decompInit();
    void decompX();
    void decompY();
    void decompZ();

    // Return processor rank from its index in the processor map
    // and vice-versa

    label procNumFromIndex(const labelVector ijk, const labelVector D);

    labelVector indexFromProcNum(const label num, const labelVector D);

    // Return processor dimensions and origins

    List<labelVector> procDims(labelVector D);
    List<labelVector> procOrig(labelVector D);

    // Unpack a receive buffer

    void unpack
    (
        const scalarBlock& buffer,
        const List<labelVector>& recvStart,
        const List<labelVector>& recvSize,
        scalarBlock& output,
        List<labelVector>& startIndex
    );

    // Transpose function

    void transpose(scalarBlock& src, scalarBlock& dst, labelVector I, labelVector T);

    // Transform Y-pencil data from i>j>k indexing to i>k>j indexing
    // and vice-versa. This is needed to apply the transformation in y-direction

    void yTransFwd(scalarBlock& yData);
    void yTransBwd(scalarBlock& yData);

};

} // end namespace briscola

} // end namespace Foam

#endif
