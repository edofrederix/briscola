#include "pointInterpolator.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
template<class Type>
List<Type> pointInterpolator<MeshType>::combine(List<Type>& values) const
{
    if (!Pstream::parRun())
        return values;

    if (global_)
    {
        // Collect at master

        DynamicList<Type> myValues;

        forAll(values, i)
            if (providers_[i] == Pstream::myProcNo())
                myValues.append(values[i]);

        List<List<Type>> all(Pstream::nProcs());
        all[Pstream::myProcNo()] = myValues;

        Pstream::gatherList(all);

        if (Pstream::master())
        {
            // Collect all processor data into the value list

            labelList c(Pstream::nProcs(), 0);

            forAll(points_, i)
                if (providers_[i] > -1)
                    values[i] = all[providers_[i]][c[providers_[i]]++];
        }

        // Distribute

        Pstream::scatter(values);
    }
    else
    {
        List<DynamicList<Type>> data(Pstream::nProcs());
        labelList recvCounts(Pstream::nProcs(), 0);

        const label me = Pstream::myProcNo();

        forAll(points_, i)
        {
            if (providers_[i] == me && requesters_[i] != me)
                data[requesters_[i]].append(values[i]);

            if (requesters_[i] == me && providers_[i] != me)
                recvCounts[providers_[i]]++;
        }

        // Setup buffers

        List<List<Type>> buffers(Pstream::nProcs());
        forAll(buffers, proc)
            if (recvCounts[proc])
                buffers[proc].resize(recvCounts[proc]);

        // Setup parallel reads

        const label nReq = Pstream::nRequests();

        forAll(buffers, proc)
            if (buffers[proc].size())
                UIPstream::read
                (
                    Pstream::commsTypes::nonBlocking,
                    proc,
                    reinterpret_cast<char*>(buffers[proc].begin()),
                    buffers[proc].byteSize()
                );

        // Setup parallel writes

        forAll(data, proc)
            if (data[proc].size())
                UOPstream::write
                (
                    Pstream::commsTypes::nonBlocking,
                    proc,
                    reinterpret_cast<char*>(data[proc].begin()),
                    data[proc].byteSize()
                );

        UPstream::waitRequests(nReq);

        // Collect all data into the value list

        labelList c(Pstream::nProcs(), 0);

        values.resize(size_);

        forAll(values, i)
            if (providers_[i] > -1 && providers_[i] != me)
                values[i] = buffers[providers_[i]][c[providers_[i]]++];

    }

    return values;
}

}

}

}
