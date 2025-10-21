#include "IO.H"
#include "OSspecific.H"
#include "Time.H"
#include "SortableList.H"
#include "colocatedFields.H"
#include "staggeredFields.H"
#include "floatVector.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(IO, 0);

word IO::dataFileNamePrefix("briscola");
word IO::vtkHeader("# vtk DataFile Version 3.0");

template<class MeshType>
void IO::writeData(const word timeName, const label l)
{
    // Prepare points

    const meshLevel<vector,MeshType>& ccl =
        fvMsh_.template metrics<MeshType>().cellCenters()[l];

    const meshLevel<vertexVector,MeshType>& vcl =
        fvMsh_.template metrics<MeshType>().vertexCenters()[l];

    const bool ascii = (fvMsh_.time().writeFormat() == IOstream::ASCII);

    forAll(ccl, d)
    {
        List<autoPtr<std::ofstream>> filePtrs;

        if (partitioned_)
        {
            filePtrs.setSize(Pstream::nProcs());
        }
        else
        {
            filePtrs.setSize(1);
        }

        forAll(filePtrs, proc)
        {
            const word vtkFile = vtkFileName<MeshType>(l,d,proc);

            const fileName vtkFilePath
            (
                fvMsh_.time().path()/timeName/vtkFile
            );

            // Create file stream on master

            if (Pstream::master())
                filePtrs[proc].reset(new std::ofstream(vtkFilePath));

            // Write header

            if (Pstream::master())
            {
                std::ofstream& file = filePtrs[proc]();

                file<< vtkHeader << std::endl
                    << fvMsh_.time().caseName() << " "
                    << (partitioned_ ? 1 : Pstream::nProcs()) << std::endl;

                if (ascii)
                {
                    file<< "ASCII" << std::endl;
                }
                else
                {
                    file<< "BINARY" << std::endl;
                }

                file<< "DATASET UNSTRUCTURED_GRID" << std::endl;
            }
        }

        // Create points data

        const labelVector S(ccl[d].I().lower()-unitXYZ*label(ghosts_));
        const labelVector E(ccl[d].I().upper()+unitXYZ*label(ghosts_));

        List<floatScalar> myPoints(cmptProduct(E-S+unitXYZ)*3);

        label c = 0;
        for (int i = S.x(); i < E.x() + 1; i++)
        for (int j = S.y(); j < E.y() + 1; j++)
        for (int k = S.z(); k < E.z() + 1; k++)
        {
            const labelVector ijk
            (
                i - (i == E.x()),
                j - (j == E.y()),
                k - (k == E.z())
            );

            const label v =
                (k == E.z())*4
              + (j == E.y())*2
              + (i == E.x());

            const floatVector p = floatVector(vcl[d](ijk)[v]);

            myPoints[c++] = p.x();
            myPoints[c++] = p.y();
            myPoints[c++] = p.z();
        }

        // Write points data

        List<labelVector> shapes(Pstream::nProcs());
        shapes[Pstream::myProcNo()] = E-S;
        Pstream::gatherList(shapes);

        forAll(filePtrs, proc)
        {
            if (Pstream::master())
            {
                label nPoints = 0;

                forAll(shapes, p)
                    if (!partitioned_ || p == proc)
                        nPoints += cmptProduct(shapes[p]+unitXYZ);

                std::ofstream& file = filePtrs[proc]();

                file<< "POINTS " << nPoints << " float" << std::endl;
            }
        }

        IO::writeList
        (
            filePtrs,
            myPoints,
            word("points"),
            1,
            fvMsh_.time().writeFormat() == IOstream::ASCII,
            l*MeshType::numberOfDirections + d,
            false
        );

        // Write cell data

        forAll(filePtrs, proc)
        {
            if (Pstream::master())
            {
                label nCells = 0;

                forAll(shapes, p)
                {
                    if (!partitioned_ || p == proc)
                    {
                        nCells += nStructured(shapes[p]);
                    }
                }

                std::ofstream& file = filePtrs[proc]();

                file<< "CELLS " << nCells << " " << 9*nCells << std::endl;

                label cursor = 0;

                forAll(shapes, p)
                if (!partitioned_ || p == proc)
                {
                    const label ll = shapes[p].x();
                    const label mm = shapes[p].y();
                    const label nn = shapes[p].z();

                    const label L = ll + 1;
                    const label M = mm + 1;
                    const label N = nn + 1;

                    labelList buffer(nStructured(shapes[p])*9);

                    label q = 0;
                    for(int i = 0; i < ll; i++)
                    for(int j = 0; j < mm; j++)
                    for(int k = 0; k < nn; k++)
                    if (structured(i,j,k,shapes[p]))
                    {
                        buffer[q++] = 8;
                        buffer[q++] = (cursor + (i  )*M*N + (j  )*N+k);
                        buffer[q++] = (cursor + (i+1)*M*N + (j  )*N+k);
                        buffer[q++] = (cursor + (i+1)*M*N + (j+1)*N+k);
                        buffer[q++] = (cursor + (i  )*M*N + (j+1)*N+k);
                        buffer[q++] = (cursor + (i  )*M*N + (j  )*N+k+1);
                        buffer[q++] = (cursor + (i+1)*M*N + (j  )*N+k+1);
                        buffer[q++] = (cursor + (i+1)*M*N + (j+1)*N+k+1);
                        buffer[q++] = (cursor + (i  )*M*N + (j+1)*N+k+1);
                    }

                    if (ascii)
                    {
                        forAll(buffer, i)
                        {
                            file<< buffer[i] << " ";
                        }
                    }
                    else
                    {
                        #ifdef LITTLEENDIAN
                        swapWords
                        (
                            buffer.size(),
                            reinterpret_cast<label*>(buffer.begin())
                        );
                        #endif

                        file.write
                        (
                            reinterpret_cast<char*>(buffer.begin()),
                            buffer.byteSize()
                        );
                    }

                    cursor += L*M*N;
                }

                file<< std::endl;

                file<< "CELL_TYPES " << nCells << std::endl;

                if (ascii)
                {
                    for (int i = 0; i < nCells; i++)
                    {
                        file<< "12 ";
                    }
                }
                else
                {
                    labelList buffer(nCells, 12);

                    #ifdef LITTLEENDIAN
                    swapWords
                    (
                        buffer.size(),
                        reinterpret_cast<label*>(buffer.begin())
                    );
                    #endif

                    file.write
                    (
                        reinterpret_cast<char*>(buffer.begin()),
                        buffer.byteSize()
                    );
                }

                file<< std::endl;
                file<< "CELL_DATA " << nCells << std::endl;
            }
        }

        // Write fields as cell data. Order is important.

        writeFields<label,MeshType>(filePtrs, l, d);
        writeFields<scalar,MeshType>(filePtrs, l, d);

        writeFields<vector,MeshType>(filePtrs, l, d);
        writeFields<tensor,MeshType>(filePtrs, l, d);
        writeFields<diagTensor,MeshType>(filePtrs, l, d);
        writeFields<sphericalTensor,MeshType>(filePtrs, l, d);
        writeFields<symmTensor,MeshType>(filePtrs, l, d);

        writeFields<faceScalar,MeshType>(filePtrs, l, d);
        writeFields<edgeScalar,MeshType>(filePtrs, l, d);
        writeFields<vertexScalar,MeshType>(filePtrs, l, d);

        writeFields<faceVector,MeshType>(filePtrs, l, d);
        writeFields<edgeVector,MeshType>(filePtrs, l, d);
        writeFields<vertexVector,MeshType>(filePtrs, l, d);
    }
}

template<class MeshType>
void IO::readData(const word timeName, const label l)
{
    const meshLevel<vector,MeshType>& ccl =
        fvMsh_.template metrics<MeshType>().cellCenters()[l];

    wordList fieldNames;

    for (label d = 0; d < MeshType::numberOfDirections; d++)
    {
        const word vtkFile = vtkFileName<MeshType>(l,d,Pstream::myProcNo());

        const fileName vtkFilePath
        (
            fvMsh_.time().path()/timeName/vtkFile
        );

        // Create file stream

        autoPtr<std::ifstream> filePtr;

        filePtr.reset(new std::ifstream(vtkFilePath));

        std::ifstream& file = filePtr();

        if (file.fail())
        {
            return;
        }

        // Read header

        char buffer[256];

        file.getline(buffer, sizeof(buffer));

        word header(buffer);

        if (header != vtkHeader && Pstream::master())
        {
            FatalErrorInFunction
                << "Invalid header in file " << vtkFilePath
                << endl << abort(FatalError);
        }

        // Case name and number of procs

        word dummy;
        label nProcs;

        file>> dummy >> nProcs;
        nextLine(file);

        if
        (
            nProcs != (partitioned_ ? 1 : Pstream::nProcs())
         && Pstream::master()
        )
            FatalErrorInFunction
                << "The vtk file " << vtkFilePath << nl
                << "was generated by " << nProcs << " processors "
                << "which is not compatible with the current number "
                << "of processors (" << Pstream::nProcs() << ")"
                << endl << abort(FatalError);

        // Ascii or binary

        word dataType;

        file>> dataType;
        nextLine(file);

        bool ascii = false;

        if (dataType == "ASCII")
        {
            ascii = true;
        }
        else if (dataType == "BINARY")
        {
            ascii = false;
        }
        else if (Pstream::master())
        {
            FatalErrorInFunction
                << "Could not determine if file " << vtkFilePath << nl
                << "is binary or ascii"
                << endl << abort(FatalError);
        }

        // Skip dataset

        nextLine(file);

        // Skip points

        if (ascii)
        {
            nextLine(file);
            nextLine(file);
        }
        else
        {
            label nPoints;

            file>> dummy >> nPoints;

            nextLine(file);
            file.ignore(nPoints*3*sizeof(floatScalar));
            nextLine(file);
        }

        // Skip cells

        if (ascii)
        {
            nextLine(file);
            nextLine(file);
        }
        else
        {
            label n;

            file>> dummy >> dummy >> n;

            nextLine(file);
            file.ignore(n*sizeof(label));
            nextLine(file);
        }

        // Skip cell types

        if (ascii)
        {
            nextLine(file);
            nextLine(file);
        }
        else
        {
            label nCells;

            file>> dummy >> nCells;

            nextLine(file);
            file.ignore(nCells*sizeof(label));
            nextLine(file);
        }

        // Read number of cells

        label nCells;

        file>> dummy >> nCells;

        const labelVector S(ccl[d].I().lower()-unitXYZ*label(ghosts_));
        const labelVector E(ccl[d].I().upper()+unitXYZ*label(ghosts_));

        labelList sizes(Pstream::nProcs());
        sizes[Pstream::myProcNo()] = nStructured(E-S);
        Pstream::gatherList(sizes);
        Pstream::scatterList(sizes);

        const label size =
            partitioned_ ? sizes[Pstream::myProcNo()] : sum(sizes);

        if (size != nCells && Pstream::master())
            FatalErrorInFunction
                << "Trying to read " << size << " cells "
                << "but found data for " << nCells << " cells in " << nl
                << vtkFilePath << endl << abort(FatalError);

        // Read fields from cell data. Order is important.

        wordList fieldNames2;

        readFields<label,MeshType>(filePtr, ascii, l, d, fieldNames2);
        readFields<scalar,MeshType>(filePtr, ascii, l, d, fieldNames2);

        readFields<vector,MeshType>(filePtr, ascii, l, d, fieldNames2);
        readFields<tensor,MeshType>(filePtr, ascii, l, d, fieldNames2);
        readFields<diagTensor,MeshType>(filePtr, ascii, l, d, fieldNames2);
        readFields<sphericalTensor,MeshType>(filePtr, ascii, l, d, fieldNames2);
        readFields<symmTensor,MeshType>(filePtr, ascii, l, d, fieldNames2);

        readFields<faceScalar,MeshType>(filePtr, ascii, l, d, fieldNames2);
        readFields<edgeScalar,MeshType>(filePtr, ascii, l, d, fieldNames2);
        readFields<vertexScalar,MeshType>(filePtr, ascii, l, d, fieldNames2);

        readFields<faceVector,MeshType>(filePtr, ascii, l, d, fieldNames2);
        readFields<edgeVector,MeshType>(filePtr, ascii, l, d, fieldNames2);
        readFields<vertexVector,MeshType>(filePtr, ascii, l, d, fieldNames2);

        if (d == 0)
            fieldNames = fieldNames2;
    }

    if (!ghosts_)
        this->correctBoundaryConditions<MeshType>(fieldNames, l);
}

template<class Type, class MeshType>
void IO::writeFields
(
    List<autoPtr<std::ofstream>>& filePtrs,
    const label l,
    const label d
) const
{
    typedef meshField<Type,MeshType> FieldType;

    const wordList toc(fvMsh_.db().toc<FieldType>());

    label count = 0;

    forAll(toc, i)
    {
        const FieldType& field =
            fvMsh_.db().lookupObject<FieldType>(toc[i]);

        if
        (
            field.writeOpt() == IOobject::AUTO_WRITE
         && (l == 0 || field.deep())
        )
        {
            count++;
        }
    }

    forAll(filePtrs, proc)
    {
        if (Pstream::master())
        {
            std::ofstream& file = filePtrs[proc]();

            file<< "FIELD " << FieldType::typeName << "s "
                << count << std::endl;
        }
    }

    if (count > 0)
    {
        forAll(toc, i)
        {
            const FieldType& field =
                fvMsh_.db().lookupObject<FieldType>(toc[i]);

            if
            (
                field.writeOpt() == IOobject::AUTO_WRITE
             && (l == 0 || field.deep())
            )
            {
                writeField(filePtrs, field[l][d]);
            }
        }
    }
}

template<class Type, class MeshType>
void IO::readFields
(
    autoPtr<std::ifstream>& filePtr,
    const bool ascii,
    const label l,
    const label d,
    wordList& fields
) const
{
    typedef meshField<Type, MeshType> FieldType;

    word dummy, type;
    label nf;

    std::ifstream& file = filePtr();

    file>> dummy >> type >> nf;

    if (word(FieldType::typeName)+"s" != type && Pstream::master())
    {
        FatalErrorInFunction
            << type << endl
            << "Could not find data for field type " << FieldType::typeName
            << endl << abort(FatalError);
    }

    wordList typeFields;

    if (nf > 0)
    {
        word name;

        for (int i = 0; i < nf; i++)
        {
            label nComponents;
            label nCells;

            file>> name >> nComponents >> nCells;
            nextLine(file);

            if (fvMsh_.db().foundObject<FieldType>(name))
            {
                FieldType& field =
                    fvMsh_.db().lookupObjectRef<FieldType>(name);

                if
                (
                    (
                        field.readOpt() == IOobject::MUST_READ
                     || field.readOpt() == IOobject::READ_IF_PRESENT
                    )
                 && (l == 0 || field.deep())
                )
                {
                    readField(filePtr, ascii, field[l][d]);
                    typeFields.append(field.name());
                }
                else
                {
                    if (ascii)
                    {
                        nextLine(file);
                    }
                    else
                    {
                        file.ignore(nCells*nComponents*sizeof(floatScalar));
                    }
                }
            }
            else
            {
                // Ignore. Alternatively we could create a field and store in
                // the registry.

                if (ascii)
                {
                    nextLine(file);
                }
                else
                {
                    file.ignore(nCells*nComponents*sizeof(floatScalar));
                }
            }
        }
    }

    HashTable<const FieldType*> objects(fvMsh_.db().lookupClass<FieldType>());

    if (l == 0)
    {
        forAllIter
        (
            typename HashTable<const FieldType*>,
            objects,
            iter
        )
        {
            if (iter()->readOpt() == IOobject::MUST_READ)
                if (findIndex(typeFields, iter()->name()) < 0)
                    FatalErrorInFunction
                        << "Count not read field " << iter()->name()
                        << " of type " << FieldType::typeName
                        << endl << abort(FatalError);
        }
    }

    fields.append(typeFields);
}

template<class Type, class MeshType>
void IO::correctBoundaryConditions(const word& fieldName, const label l)
{
    if (fvMsh_.db().foundObject<meshField<Type,MeshType>>(fieldName))
    {
        auto& field =
            fvMsh_.db().lookupObjectRef<meshField<Type,MeshType>>(fieldName);

        field[l].correctBoundaryConditions();
    }
}

template<class MeshType>
void IO::correctBoundaryConditions(const wordList& fieldNames, const label l)
{
    forAll(fieldNames, i)
    {
        const word fieldName = fieldNames[i];

        correctBoundaryConditions<label,MeshType>(fieldName, l);
        correctBoundaryConditions<scalar,MeshType>(fieldName, l);
        correctBoundaryConditions<vector,MeshType>(fieldName, l);
        correctBoundaryConditions<tensor,MeshType>(fieldName, l);
        correctBoundaryConditions<diagTensor,MeshType>(fieldName, l);
        correctBoundaryConditions<sphericalTensor,MeshType>(fieldName, l);
        correctBoundaryConditions<symmTensor,MeshType>(fieldName, l);
        correctBoundaryConditions<faceScalar,MeshType>(fieldName, l);
        correctBoundaryConditions<edgeScalar,MeshType>(fieldName, l);
        correctBoundaryConditions<vertexScalar,MeshType>(fieldName, l);
        correctBoundaryConditions<faceVector,MeshType>(fieldName, l);
        correctBoundaryConditions<edgeVector,MeshType>(fieldName, l);
        correctBoundaryConditions<vertexVector,MeshType>(fieldName, l);
    }
}

template<class MeshType>
word IO::vtkFileName
(
    const label l,
    const label d,
    const label p
) const
{
    return
        dataFileNamePrefix
      + "_"
      + (l > 0 ? word(Foam::name(l) + "_") : "")
      + (partitioned_ ? word(Foam::name(p) + "_") : "")
      + MeshType::typeName
      + (
            MeshType::numberOfDirections > 1
          ? word("_" + word(MeshType::directionNames[d]))
          : word("")
        )
      + word(".vtk");
}

template<class MeshType>
word IO::seriesFileName(const label l, const label d, const label p) const
{
    return
        dataFileNamePrefix
      + "_"
      + (l > 0 ? word(Foam::name(l) + "_") : "")
      + (partitioned_ ? word(Foam::name(p) + "_") : "")
      + MeshType::typeName
      + (
            MeshType::numberOfDirections > 1
          ? word("_" + word(MeshType::directionNames[d]))
          : word("")
        )
      + ".vtk.series";
}

template<class MeshType>
void IO::writeSeries(const label l) const
{
    if (Pstream::master())
    {
        instantList times(fvMsh_.time().times());

        for (label p = 0; p < (partitioned_ ? Pstream::nProcs() : 1); p++)
        for (label d = 0; d < MeshType::numberOfDirections; d++)
        {
            OFstream seriesFile
            (
                fvMsh_.time().path()/seriesFileName<MeshType>(l,d,p)
            );

            seriesFile
                << "{" << nl
                << "    \"file-series-version\" : \"1.0\"," << nl
                << "    \"files\" :" << nl
                << "    [" << nl;

            forAll(times, i)
            if (times[i].value() != 0)
            {
                seriesFile
                    << "        { \"name\" : \""
                    << times[i].name() << "/" << vtkFileName<MeshType>(l,d,p)
                    << "\", \"time\" : "
                    << times[i].name() << " }," << nl;
            }

            seriesFile
                << "    ]" << nl
                << "}" << endl;
        }
    }
}

void IO::swapWord(label& word32) const
{
    char* mem = reinterpret_cast<char*>(&word32);

    char a = mem[0];
    mem[0] = mem[3];
    mem[3] = a;

    a = mem[1];
    mem[1] = mem[2];
    mem[2] = a;
}


void IO::swapWords(const label nWords, label* words32) const
{
    for (label i = 0; i < nWords; i++)
    {
        swapWord(words32[i]);
    }
}

IO::IO(const fvMesh& fvMsh)
:
    fvMsh_(fvMsh),
    partitioned_
    (
        fvMsh.time().controlDict()
       .lookupOrDefault<Switch>("writePartitioned", false)
    ),
    ghosts_
    (
        fvMsh.time().controlDict()
       .lookupOrDefault<Switch>("writeGhosts", false)
    )
{}

IO::~IO()
{}

template<class MeshType>
void IO::write(const label l)
{
    if (fvMsh_.time().writeTime())
    {
        writeNow<MeshType>(l);
    }
}

template<class MeshType>
void IO::writeNow(const label l)
{
    const Time& time = fvMsh_.time();

    const word timeName(time.name());

    if (Pstream::master())
    {
        if (!isDir(time.path()/timeName))
        {
            mkDir(time.path()/timeName);
        }
    }

    // Wait

    if (Pstream::parRun())
        returnReduce(1.0,sumOp<scalar>());

    // Write data

    writeData<MeshType>(timeName,l);

    // Write series

    writeSeries<MeshType>(l);
}

template<class MeshType>
void IO::read(const label l)
{
    const Time& time = fvMsh_.time();

    const word timeName(time.name());

    if (Pstream::master())
    {
        if (!isDir(time.path()/timeName) || timeName == "0")
        {
            FatalErrorInFunction
                << "Cannot read from time " << timeName << endl
                << exit(FatalError);
        }
    }

    // Read data

    readData<MeshType>(timeName,l);
}

// Instantiate

template void IO::read<colocated>(const label);
template void IO::read<staggered>(const label);

template void IO::write<colocated>(const label);
template void IO::write<staggered>(const label);

}

}

}
