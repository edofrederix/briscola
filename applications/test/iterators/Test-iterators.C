#include "arguments.H"
#include "Time.H"

#include "fvMesh.H"
#include "faceFields.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

// Mesh is supposed to be 8*16*32

const labelVector shape(8,16,32);

template<class Type, class MeshType>
void testCellIterators(const meshField<Type,MeshType>& field)
{
    label refCount = 0;
    forAll(field[0], d)
        refCount += cmptProduct(shape + MeshType::padding[d]);

    label count;

    #define CELLCHECK(d,i,j,k)                                                 \
        i < 0 || i >= (shape.x() + MeshType::padding[d].x())                   \
     || j < 0 || j >= (shape.y() + MeshType::padding[d].y())                   \
     || k < 0 || k >= (shape.z() + MeshType::padding[d].z())

    // Direction

    count = 0;

    forAll(field[0], d)
    forAllCellsInDirection(field[0][d], i, j, k)
    {
        if (CELLCHECK(d,i,j,k))
            FatalErrorInFunction
                << "Invalid cell index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not all cells were traversed" << endl << abort(FatalError);

    count = 0;

    forAll(field[0], d)
    forAllCells(field[0][d], i, j, k)
    {
        if (CELLCHECK(d,i,j,k))
            FatalErrorInFunction
                << "Invalid cell index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not all cells were traversed" << endl << abort(FatalError);

    // Level

    count = 0;

    forAllCellsInLevel(field[0], d, i, j, k)
    {
        if (CELLCHECK(d,i,j,k))
            FatalErrorInFunction
                << "Invalid cell index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not all cells were traversed" << endl << abort(FatalError);

    count = 0;

    forAllCells(field[0], d, i, j, k)
    {
        if (CELLCHECK(d,i,j,k))
            FatalErrorInFunction
                << "Invalid cell index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not all cells were traversed" << endl << abort(FatalError);

    // Field

    count = 0;

    forAllCellsInField(field, l, d, i, j, k)
    {
        if (CELLCHECK(d,i,j,k))
            FatalErrorInFunction
                << "Invalid cell index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not all cells were traversed" << endl << abort(FatalError);

    count = 0;

    forAllCells(field, l, d, i, j, k)
    {
        if (CELLCHECK(d,i,j,k))
            FatalErrorInFunction
                << "Invalid cell index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not all cells were traversed" << endl << abort(FatalError);

    // Reversed direction

    count = 0;

    forAll(field[0], d)
    forAllCellsReversedInDirection(field[0][d], i, j, k)
    {
        if (CELLCHECK(d,i,j,k))
            FatalErrorInFunction
                << "Invalid cell index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not all cells were traversed" << endl << abort(FatalError);

    count = 0;

    forAll(field[0], d)
    forAllCellsReversed(field[0][d], i, j, k)
    {
        if (CELLCHECK(d,i,j,k))
            FatalErrorInFunction
                << "Invalid cell index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not all cells were traversed" << endl << abort(FatalError);

    // Reversed level

    count = 0;

    forAllCellsReversedInLevel(field[0], d, i, j, k)
    {
        if (CELLCHECK(d,i,j,k))
            FatalErrorInFunction
                << "Invalid cell index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not all cells were traversed" << endl << abort(FatalError);

    count = 0;

    forAllCellsReversed(field[0], d, i, j, k)
    {
        if (CELLCHECK(d,i,j,k))
            FatalErrorInFunction
                << "Invalid cell index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not all cells were traversed" << endl << abort(FatalError);

    // Reversed field

    count = 0;

    forAllCellsReversedInField(field, l, d, i, j, k)
    {
        if (CELLCHECK(d,i,j,k))
            FatalErrorInFunction
                << "Invalid cell index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not all cells were traversed" << endl << abort(FatalError);

    count = 0;

    forAllCellsReversed(field, l, d, i, j, k)
    {
        if (CELLCHECK(d,i,j,k))
            FatalErrorInFunction
                << "Invalid cell index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not all cells were traversed" << endl << abort(FatalError);
}

template<class Type, class MeshType>
void testFaceIterators(const faceField<Type,MeshType>& set)
{
    using briscola::units;

    label refCount = 0;
    forAll(set, fd)
        forAll(set[fd][0], d)
            refCount +=
                cmptProduct(shape + MeshType::padding[d] + units[fd]);

    label count;

    #define FACECHECK(fd,d,i,j,k)                                              \
        i < 0 || i >= (shape.x() + MeshType::padding[d].x() + units[fd].x())   \
     || j < 0 || j >= (shape.y() + MeshType::padding[d].y() + units[fd].y())   \
     || k < 0 || k >= (shape.z() + MeshType::padding[d].z() + units[fd].z())

     // Direction

    count = 0;

    forAll(set[0][0], d)
    forAllFacesInSpecificDirection(set, fd, 0, d, i, j, k)
    {
        if (FACECHECK(fd,d,i,j,k))
            FatalErrorInFunction
                << "Invalid face index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not all faces were traversed" << endl << abort(FatalError);

    // Level

    count = 0;

    forAllFacesInLevel(set, fd, d, i, j, k)
    {
        if (FACECHECK(fd,d,i,j,k))
            FatalErrorInFunction
                << "Invalid face index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not all cells were traversed" << endl << abort(FatalError);

    // Field

    count = 0;

    forAllFacesInField(set, fd, l, d, i, j, k)
    {
        if (FACECHECK(fd,d,i,j,k))
            FatalErrorInFunction
                << "Invalid face index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not all faces were traversed" << endl << abort(FatalError);

    // General

    count = 0;

    forAllFaces(set, fd, l, d, i, j, k)
    {
        if (FACECHECK(fd,d,i,j,k))
            FatalErrorInFunction
                << "Invalid face index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not face cells were traversed" << endl << abort(FatalError);
}

template<class Type, class MeshType>
void testInternalFaceIterators(const faceField<Type,MeshType>& set)
{
    using briscola::units;

    label refCount = 0;
    forAll(set, fd)
        forAll(set[fd][0], d)
            refCount +=
                cmptProduct(shape + MeshType::padding[d] - units[fd]);

    label count;

    #define INTERNALFACECHECK(fd,d,i,j,k)                                      \
        i < units[fd].x()                                                      \
     || j < units[fd].y()                                                      \
     || k < units[fd].z()                                                      \
     || i >= (shape.x() + MeshType::padding[d].x())                            \
     || j >= (shape.y() + MeshType::padding[d].y())                            \
     || k >= (shape.z() + MeshType::padding[d].z())

    // Level

    count = 0;

    forAllInternalFacesInLevel(set, fd, d, i, j, k)
    {
        if (INTERNALFACECHECK(fd,d,i,j,k))
            FatalErrorInFunction
                << "Invalid face index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not all cells were traversed" << endl << abort(FatalError);

    // Field

    count = 0;

    forAllInternalFacesInField(set, fd, l, d, i, j, k)
    {
        if (INTERNALFACECHECK(fd,d,i,j,k))
            FatalErrorInFunction
                << "Invalid face index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not all faces were traversed" << endl << abort(FatalError);

    count = 0;

    forAllInternalFaces(set, fd, l, d, i, j, k)
    {
        if (INTERNALFACECHECK(fd,d,i,j,k))
            FatalErrorInFunction
                << "Invalid face index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not face cells were traversed" << endl << abort(FatalError);
}

template<class Type, class MeshType>
void testBoundaryFaceIterators(const faceField<Type,MeshType>& set)
{
    using briscola::units;

    label refCount = 0;
    forAll(set, fd)
        forAll(set[fd][0], d)
            refCount +=
                2
              * cmptProduct(shape + MeshType::padding[d])
              / (shape + MeshType::padding[d])[fd];

    label count;

    #define BOUNDARYFACECHECK(fd,d,i,j,k)                                      \
        (                                                                      \
            fd == 0                                                            \
         && (                                                                  \
                (i != 0 && i != shape.x() + MeshType::padding[d].x())          \
             || (j < 0  || j >= shape.y() + MeshType::padding[d].y())          \
             || (k < 0  || k >= shape.z() + MeshType::padding[d].z())          \
            )                                                                  \
        )                                                                      \
     || (                                                                      \
            fd == 1                                                            \
         && (                                                                  \
                (j != 0 && j != shape.y() + MeshType::padding[d].y())          \
             || (i < 0  || i >= shape.x() + MeshType::padding[d].x())          \
             || (k < 0  || k >= shape.z() + MeshType::padding[d].z())          \
            )                                                                  \
        )                                                                      \
     || (                                                                      \
            fd == 2                                                            \
         && (                                                                  \
                (k != 0 && k != shape.z() + MeshType::padding[d].z())          \
             || (i < 0  || i >= shape.x() + MeshType::padding[d].x())          \
             || (j < 0  || j >= shape.y() + MeshType::padding[d].y())          \
            )                                                                  \
        )

    // Level

    count = 0;

    forAllBoundaryFacesInLevel(set, fd, lu, d, i, j, k)
    {
        if (BOUNDARYFACECHECK(fd,d,i,j,k))
            FatalErrorInFunction
                << "Invalid face index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not all cells were traversed" << endl << abort(FatalError);

    // Field

    count = 0;

    forAllBoundaryFacesInField(set, fd, lu, l, d, i, j, k)
    {
        if (BOUNDARYFACECHECK(fd,d,i,j,k))
            FatalErrorInFunction
                << "Invalid face index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not all faces were traversed" << endl << abort(FatalError);

    count = 0;

    forAllBoundaryFaces(set, fd, lu, l, d, i, j, k)
    {
        if (BOUNDARYFACECHECK(fd,d,i,j,k))
            FatalErrorInFunction
                << "Invalid face index" << endl << abort(FatalError);

        count++;
    }

    if (count != refCount)
        FatalErrorInFunction
            << "Not face cells were traversed" << endl << abort(FatalError);
}

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"

    const colocatedScalarField coloCells
    (
        colocatedScalarField::New("field", fvMsh)
    );

    const staggeredScalarField stagCells
    (
        staggeredScalarField::New("field", fvMsh)
    );

    const colocatedScalarFaceField coloFaces
    (
        colocatedScalarFaceField::New("field", fvMsh)
    );

    const staggeredScalarFaceField stagFaces
    (
        staggeredScalarFaceField::New("field", fvMsh)
    );

    testCellIterators(coloCells);
    testCellIterators(stagCells);

    testFaceIterators(coloFaces);
    testFaceIterators(stagFaces);

    testInternalFaceIterators(coloFaces);
    testInternalFaceIterators(stagFaces);

    testBoundaryFaceIterators(coloFaces);
    testBoundaryFaceIterators(stagFaces);
}
