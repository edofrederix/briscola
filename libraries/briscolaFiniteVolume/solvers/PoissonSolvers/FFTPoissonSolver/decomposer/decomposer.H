#ifndef decomposer_H
#define decomposer_H

#include "vector.H"
#include "block.H"
#include "word.H"
#include "fvMesh.H"
#include "rectilinearMesh.H"
#include "processorOverlap.H"
#include "HashTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class decomposer
{
	// Mesh
    const fvMesh& fvMsh_;

    // Global Size
	labelVector N_;

    // Initial decomposition
    labelVector I_;

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
    List<labelVector> Si_;
    List<labelVector> Sx_;
    List<labelVector> Sy_;
    List<labelVector> Sz_;

    // Processor overlap lookup table
    HashTable<processorOverlap> processorOverlap_;

public:

    // Construct from mesh
    decomposer(const fvMesh& fvMsh, label decompType);

    // Destructor
    virtual ~decomposer();

    // Return initial decomposition
    inline labelVector& I() { return I_; }
    inline List<labelVector>& Ni() { return Ni_; }
    inline List<labelVector>& Si() { return Si_; }

    // Return X-pencil decomposition
    inline labelVector& X() { return X_; }
    inline List<labelVector>& Nx() { return Nx_; }
    inline List<labelVector>& Sx() { return Sx_; }

    // Return Y-pencil decomposition
    inline labelVector& Y() { return Y_; }
    inline List<labelVector>& Ny() { return Ny_; }
    inline List<labelVector>& Sy() { return Sy_; }

    // Return Z-pencil decomposition
    inline labelVector& Z() { return Z_; }
    inline List<labelVector>& Nz() { return Nz_; }
    inline List<labelVector>& Sz() { return Sz_; }

    // Return pencil decomposition for a given direction
    List<labelVector>& Nd(label dir)
    {
        switch (dir)
        {
            case 0:
                return Nx_;
                break;

            case 1:
                return Ny_;
                break;

            case 2:
                return Nz_;
                break;

            default:
                FatalError
                    << "Invalid direction."
                    << endl;
                FatalError.exit();
                return Nx_;
                break;
        }
    }

    List<labelVector>& Sd(label dir)
    {
        switch (dir)
        {
            case 0:
                return Sx_;
                break;

            case 1:
                return Sy_;
                break;

            case 2:
                return Sz_;
                break;

            default:
                FatalError
                    << "Invalid direction."
                    << endl;
                FatalError.exit();
                return Sx_;
                break;
        }
    }

    // Initialize decompositions
    void decompInit(label decompType);

    // Check a given decomposition
    void checkDecomp(labelVector D);

    // Return processor rank from its index in the processor map
    // and vice-versa
    label procNumFromIndex(const labelVector ijk, const labelVector D);
    labelVector indexFromProcNum(const label num, const labelVector D);

    // Make key for HashTable
    word makeKey(labelVector I, labelVector T);

    // Return processor dimensions and origins
    List<labelVector> procDims(labelVector D);
    List<labelVector> procOrig(labelVector D);

    // Unpack a receive buffer
    // recvMajorOrder represents the major order
    // of the output block ("x", "y" or "z")
    void unpack
    (
        const scalarBlock& buffer,
        const List<labelVector>& recvStart,
        const List<labelVector>& recvSize,
        scalarBlock& output,
        List<labelVector>& startIndex,
        labelVector T,
        string recvMajorOrder
    );

    // Transpose function
    // sendMajorOrder represents the major order
    // of the source block ("x", "y" or "z")
    void transpose
    (
        const scalarBlock& src,
        scalarBlock& dst,
        labelVector I,
        labelVector T,
        string sendMajorOrder,
        string recvMajorOrder
    );
};

} // end namespace fv

} // end namespace briscola

} // end namespace Foam

#endif
