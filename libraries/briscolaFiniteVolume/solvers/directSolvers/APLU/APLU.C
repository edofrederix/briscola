#include "APLU.H"

#include "SquareMatrix.H"
#include "LUscalarMatrix.H"
#include "meshDirectionStencilFunctions.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
APLU<SType,Type,MeshType>::APLU
(
    const dictionary& dict,
    const fvMesh& fvMsh,
    const label l
)
:
    directSolver<SType,Type,MeshType>(dict,fvMsh,l),
    shapes_(MeshType::numberOfDirections),
    sizes_(MeshType::numberOfDirections),
    commSizes_(Pstream::nProcs(),0),
    indices_(MeshType::numberOfDirections),
    n_(MeshType::numberOfDirections, 0)
{
    const decomposition& decomp = fvMsh.msh().decomp();

    const List<faceLabel>& fNeigh = decomp.faceNeighborsPerProc();
    const List<edgeLabel>& eNeigh = decomp.edgeNeighborsPerProc();
    const List<vertexLabel>& vNeigh = decomp.vertexNeighborsPerProc();

    const List<faceLabelTensor>& fT = decomp.faceTsPerProc();
    const List<edgeLabelTensor>& eT = decomp.edgeTsPerProc();
    const List<vertexLabelTensor>& vT = decomp.vertexTsPerProc();

    for (int d = 0; d < MeshType::numberOfDirections; d++)
    {
        // Collect data shapes

        shapes_[d].setSize(Pstream::nProcs());

        List<labelVector>& shapes = shapes_[d];

        shapes[Pstream::myProcNo()] = fvMsh.N<MeshType>(l,d);
        Pstream::gatherList(shapes);

        commSizes_[Pstream::myProcNo()] +=
            cmptProduct(shapes[Pstream::myProcNo()]);

        // On master, determine the global indices of the full stencil elements
        // of all cells

        if (Pstream::master())
        {
            // Set sizes

            sizes_[d].setSize(Pstream::nProcs());
            labelList& sizes = sizes_[d];

            forAll(sizes, proc)
                sizes[proc] = cmptProduct(shapes[proc]);

            indices_[d].setSize(Pstream::nProcs());

            // Index displacements per processor

            labelList disp(Pstream::nProcs(),0);

            disp[0] = 0;
            for (int proc = 1; proc < Pstream::nProcs(); proc++)
                disp[proc] = disp[proc-1] + sizes[proc-1];

            List<List<labelList>>& indices = indices_[d];

            forAll(indices, proc)
            {
                const labelVector shape = shapes[proc];
                const label size = sizes[proc];

                // Initialize indices at -1

                indices[proc].setSize
                (
                    size,
                    labelList(stencil::nComponents, -1)
                );

                int l = 0;
                for (int i = 0; i < shape.x(); i++)
                for (int j = 0; j < shape.y(); j++)
                for (int k = 0; k < shape.z(); k++)
                {
                    const labelVector ijk(i,j,k);

                    for (int s = 0; s < stencil::nComponents; s++)
                    {
                        const labelVector offset(stencil::componentOffsets[s]);
                        labelVector cell(ijk + offset);

                        if
                        (
                            cmptMin(cell) >= 0
                         && cmptMin(shape-cell-unitXYZ) >= 0
                        )
                        {
                            // Processor-local cell

                            indices[proc][l][s] =
                                disp[proc]
                              + cell.x()*shape.y()*shape.z()
                              + cell.y()*shape.z()
                              + cell.z();
                        }
                        else
                        {
                            // Non-local cell. Calculate index if it is in
                            // another processor.

                            const labelVector bo =
                                briscola::cmptMax
                                (
                                    briscola::cmptMin(offset, unitXYZ),
                                  - unitXYZ
                                );

                            const label degree = cmptSum(cmptMag(bo));

                            label num, neigh;
                            labelTensor T;

                            if (degree == 1)
                            {
                                num = faceNumber(bo);
                                neigh = fNeigh[proc][num];
                                T = fT[proc][num];
                            }
                            else if (degree == 2)
                            {
                                num = edgeNumber(bo);
                                neigh = eNeigh[proc][num];
                                T = eT[proc][num];
                            }
                            else
                            {
                                num = vertexNumber(bo);
                                neigh = vNeigh[proc][num];
                                T = vT[proc][num];
                            }

                            if (neigh > -1)
                            {
                                // Get the shape of the neighboring processor
                                // and transform it to the current processor's
                                // local coordinate system

                                const labelVector shape = shapes[neigh];
                                const labelVector shapeT = cmptMag(T & shape);

                                // Determine the cell in the neighbor's
                                // processor

                                for (int f = 0; f < 3; f++)
                                    if (bo[f] != 0)
                                        cell[f] =
                                            offset[f] < 0
                                          ? offset[f] + shapeT[f]
                                          : offset[f] - 1;

                                // Transform back to the neighbor's coordinate
                                // system and determine index

                                cell = (T.T() & cell) + shape;

                                for (int f = 0; f < 3; f++)
                                    cell[f] = (cell[f] % shape[f]);

                                indices[proc][l][s] =
                                    disp[neigh]
                                  + cell.x()*shape.y()*shape.z()
                                  + cell.y()*shape.z()
                                  + cell.z();
                            }
                        }
                    }

                    l++;
                }
            }
        }
    }

    // Collect total communication sizes

    Pstream::gatherList(commSizes_);

    // Compute system sizes

    if (Pstream::master())
        forAll(shapes_, d)
            forAll(shapes_[d], proc)
                n_[d] += sizes_[d][proc];
}

template<class SType, class Type, class MeshType>
void APLU<SType,Type,MeshType>::solve
(
    linearSystem<SType,Type,MeshType>& xEqn,
    const List<bool>& singular
)
{
    // Prepare data, collecting all directions in a single list

    List<Row<stencil,Type>> myData(commSizes_[Pstream::myProcNo()]);

    const meshLevel<SType,MeshType>& A = xEqn.A()[this->l_];
    const meshLevel<Type,MeshType>& b = xEqn.b()[this->l_];

    int l = 0;
    forAllCells(A, d, i, j, k)
        myData[l++] =
            Row<stencil,Type>
            (
                fullStencil<MeshType>(A[d], labelVector(i,j,k)),
                b(d,i,j,k)
            );

    List<Type> solution(myData.size());

    if (!Pstream::master())
    {
        // Send my data to master

        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            Pstream::masterNo(),
            reinterpret_cast<char*>(myData.begin()),
            myData.byteSize(),
            0,
            UPstream::worldComm
        );

        // Get solution

        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            Pstream::masterNo(),
            reinterpret_cast<char*>(solution.begin()),
            solution.byteSize(),
            0,
            UPstream::worldComm
        );
    }
    else
    {
        // Collect data from slaves

        List<List<Row<stencil,Type>>> allData(Pstream::nProcs());

        for (int proc = 0; proc < Pstream::nProcs(); proc++)
        {
            if (proc == Pstream::myProcNo())
            {
                allData[proc] = myData;
            }
            else
            {
                allData[proc].setSize(commSizes_[proc]);

                UIPstream::read
                (
                    Pstream::commsTypes::blocking,
                    proc,
                    reinterpret_cast<char*>(allData[proc].begin()),
                    allData[proc].byteSize(),
                    0,
                    UPstream::worldComm
                );
            }
        }

        // Build and solve linear systems

        labelList directionDisp(Pstream::nProcs(),0);
        List<List<Type>> solutions(Pstream::nProcs());

        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            // Create linear system. Increase dimension by one, to allow for
            // singular system augmentation.

            scalarSquareMatrix A(n_[d]+1, Zero);
            List<Type> b(n_[d]+1);
            b[n_[d]] = Zero;

            label procOffset = 0;
            forAll(indices_[d], proc)
            {
                const List<Row<stencil,Type>>& data = allData[proc];
                const List<labelList>& indices = indices_[d][proc];

                for (int l = 0; l < sizes_[d][proc]; l++)
                {
                    const label ld = directionDisp[proc] + l;

                    b[procOffset+l] = data[ld].source();

                    // When we have a non-zero stencil element but the index is
                    // -1, it means we are on a non-eliminated non-periodic and
                    // non-parallel boundary. In that case the system cannot be
                    // solved directly, because the ghost value is unknown. For
                    // now, a homogeneous Neumann boundary condition is assumed.

                    const stencil S = data[ld].stencil();
                    for (int s = 0; s < stencil::nComponents; s++)
                        if (indices[l][s] > -1)
                            A[procOffset+l][indices[l][s]] = S[s];
                        else if (Foam::mag(S[s]) > 1e-8)
                            A[procOffset+l][procOffset+l] += S[s];
                }

                procOffset += sizes_[d][proc];
            }

            if (singular[d])
            {
                // Singular matrix, augment the system (see "Multigrid", U.
                // Trottenberg et al., 2001)

                for (int i = 0; i < n_[d]; i++)
                {
                    A(i,n_[d]) = 1.0;
                    A(n_[d],i) = 1.0;
                }
            }
            else
            {
                // Set the auxilary variable to zero

                A(n_[d],n_[d]) = 1.0;
            }

            // Solve

            LUsolve(A,b);

            // Store solution

            procOffset = 0;
            forAll(solutions, proc)
            {
                if (d == 0)
                    solutions[proc].setSize(commSizes_[proc]);

                for (int l = 0; l < sizes_[d][proc]; l++)
                {
                    const label ld = directionDisp[proc] + l;

                    solutions[proc][ld] = b[procOffset+l];
                }

                procOffset += sizes_[d][proc];
            }

            // Increase data direction displacement

            if (MeshType::numberOfDirections > 1)
                forAll(directionDisp, proc)
                    directionDisp[proc] += sizes_[d][proc];
        }

        // Send to slaves

        forAll(solutions, proc)
        {
            if (proc == Pstream::myProcNo())
            {
                solution = solutions[proc];
            }
            else
            {
                UOPstream::write
                (
                    Pstream::commsTypes::blocking,
                    proc,
                    reinterpret_cast<char*>(solutions[proc].begin()),
                    solutions[proc].byteSize(),
                    0,
                    UPstream::worldComm
                );
            }
        }
    }

    // Write the solution back to the linear system

    meshLevel<Type,MeshType>& x = xEqn.x()[this->l_];

    l = 0;
    forAllCells(x, d, i, j, k)
        x(d,i,j,k) = solution[l++];

    x.correctBoundaryConditions();
}

}

}

}
