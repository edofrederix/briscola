#include "brickTopology.H"
#include "geometry.H"
#include "patch.H"

namespace Foam
{

namespace briscola
{

void brickTopology::setFaceLinks()
{
    forAll(geo_.bricks(), bricki)
    forAll(geo_.bricks(), brickj)
    if (bricki != brickj)
    {
        const brick& b0 = geo_.bricks()[bricki];
        const brick& b1 = geo_.bricks()[brickj];

        for (label facei = 0; facei < 6; facei++)
        for (label facej = 0; facej < 6; facej++)
        {
            const face& f0 = b0.f(facei);
            const face& f1 = b1.f(facej);

            if (f0 == f1)
            {
                links_[bricki].faceLinks().set
                (
                    facei,
                    new brickFaceLink(f0, f1)
                );
            }
        }
    }
}

void brickTopology::setPeriodicFaceLinks()
{
    forAll(geo_.patchPairs(), pairi)
    {
        const patchPair& pair = geo_.patchPairs()[pairi];

        const patch& patch0 = pair.patch0();
        const patch& patch1 = pair.patch1();

        forAll(patch0.facePtrs(), facei)
        {
            const face& f0 = *patch0.facePtrs()[facei];
            const face& f1 = *patch1.facePtrs()[facei];

            const brick& b0 = f0.parentBrick();
            const brick& b1 = f1.parentBrick();

            const labelList P01 = shortestFacePath(b0.num(), b1.num());

            links_[b0.num()].faceLinks().set
            (
                f0.num(),
                new brickFaceLink(f0, f1, P01, *this, true)
            );

            links_[b1.num()].faceLinks().set
            (
                f1.num(),
                new brickFaceLink(f1, f0, reverseList(P01), *this, true)
            );
        }
    }
}

void brickTopology::setEdgeAndVertexLinks()
{
    // Search edge links

    forAll(geo_.bricks(), bricki)
    {
        const brick& b0 = geo_.bricks()[bricki];

        forAll(links_[bricki].faceLinks(), facei0)
        if (links_[bricki].faceLinks().set(facei0))
        {
            const brickFaceLink& link01 = links_[bricki].faceLinks()[facei0];

            const label brickj = link01.f1().parentBrick().num();
            const label facei1 = link01.f1().num();

            // Look for face links to a third brick that are perpendicular to
            // the first link, thus forming an edge link between the first and
            // third brick

            forAll(links_[brickj].faceLinks(), facej1)
            if
            (
                facej1/2 != facei1/2
             && links_[brickj].faceLinks().set(facej1)
            )
            {
                const brickFaceLink& link12 = links_[brickj].faceLinks()[facej1];

                const label brickk = link12.f1().parentBrick().num();
                const label facej2 = link12.f1().num();

                const brick& b2 = geo_.bricks()[brickk];

                const label facei2 =
                    faceNumber(link12.T().T() & faceOffsets[facei1]);
                const label facej0 =
                    faceNumber(link01.T() & faceOffsets[facej1]);

                const labelTensor T = (link01.T() & link12.T());

                const edge& e0 = b0.e(sharedEdgeNumber(facei0,facej0));
                const edge& e1 = b2.e(sharedEdgeNumber(facei2,facej2));

                if (!links_[bricki].edgeLinks().set(e0.num()))
                {
                    links_[bricki].edgeLinks().set
                    (
                        e0.num(),
                        new brickEdgeLink(e0, e1, T, e0 != e1)
                    );
                }

                // Look for face links to a fourth brick that are perpendicular
                // to the first two links, thus forming a vertex link between
                // the first and fourth brick

                forAll(links_[brickk].faceLinks(), facek2)
                if
                (
                    facek2/2 != facei2/2
                 && facek2/2 != facej2/2
                 && links_[brickk].faceLinks().set(facek2)
                )
                {
                    const brickFaceLink& link23 = links_[brickk].faceLinks()[facek2];

                    const label brickl = link23.f1().parentBrick().num();
                    const label facek3 = link23.f1().num();

                    const brick& b3 = geo_.bricks()[brickl];

                    const label facei3 =
                        faceNumber(link23.T().T() & faceOffsets[facei2]);
                    const label facej3 =
                        faceNumber(link23.T().T() & faceOffsets[facej2]);
                    const label facek0 =
                        faceNumber(link01.T() & link12.T() & faceOffsets[facek2]);

                    const labelTensor T = (link01.T() & link12.T() & link23.T());

                    const vertex& v0 = b0.v(sharedVertexNumber(facei0,facej0,facek0));
                    const vertex& v1 = b3.v(sharedVertexNumber(facei3,facej3,facek3));

                    if (!links_[bricki].vertexLinks().set(v0.num()))
                    {
                        links_[bricki].vertexLinks().set
                        (
                            v0.num(),
                            new brickVertexLink(v0, v1, T, v0 != v1)
                        );
                    }
                }
            }
        }
    }
}

bool brickTopology::buildTopologyMap
(
    const brick& b0,
    const labelTensor T,
    labelVector& cursor,
    labelList& found,
    labelBlock& map
) const
{
    const label bricki = b0.num();

    // Update search registry and map

    found[bricki] = 1;
    map(cursor) = bricki;

    // Look at existing face links

    for (label facei = 0; facei < 6; facei++)
    if (links_[bricki].faceLinks().set(facei))
    {
        const brickFaceLink& link = links_[bricki].faceLinks()[facei];
        const brick& b1 = link.f1().parentBrick();

        // Move cursor to new position

        cursor += (T & faceOffsets[facei]);

        // Expand map if needed

        for(label dir = 0; dir < 3; dir++)
        {
            if (cursor[dir] == -1)
            {
                map.prepend(dir,1,-1);
                cursor[dir] = 0;
            }
            else if (cursor[dir] == map.shape()[dir])
            {
                map.append(dir,1,-1);
            }
        }

        // If the neighboring brick overrides another brick at the cursor, or if
        // it has already been found on another position, then the bricks are
        // unstructured. Abort in either case.

        if
        (
            (map(cursor) != -1 && map(cursor) != b1.num())
         || (map(cursor) == -1 && found[b1.num()])
        )
        {
            return false;
        }

        // If the neighboring brick is not yet found, continue building the
        // topology map for that brick. Abort if building the topology map for
        // that brick fails.

        if (!found[b1.num()])
        {
            if (!buildTopologyMap(b1, T & link.T(), cursor, found, map))
            {
                return false;
            }
        }

        // Move cursor back to the parent brick

        cursor -= (T & faceOffsets[facei]);
    }

    return true;
}

brickTopology::brickTopology(const geometry& geo)
:
    geo_(geo),
    links_(geo.bricks().size())
{
    forAll(geo.bricks(), bricki)
    {
        links_.set(bricki, new brickLinks(geo.bricks()[bricki]));
    }

    // Set face links

    setFaceLinks();

    // Try to build the topology map. If it returns false, the bricks are
    // unstructured

    labelList found(geo.bricks().size(), 0);
    labelBlock map(unitXYZ,0);
    labelVector cursor(zeroXYZ);

    structured_ = buildTopologyMap
    (
        geo.bricks()[0],
        eye,
        cursor,
        found,
        map
    );

    // Store topology map if the bricks are structured

    if (structured_)
    {
        map_.setData(map);
    }

    // If structured, check if rectilinear

    if (structured_)
    {
        rectilinear_ = min(map) > -1;
    }
    else
    {
        rectilinear_ = false;
    }

    // Check if bricks are aligned

    aligned_ = true;

    if (geo_.bricks().size() > 1)
    {
        forAll(geo_.bricks(), bricki)
        if (aligned_)
        {
            forAll(links_[bricki].faceLinks(), facei)
            if (aligned_)
            {
                if
                (
                    links_[bricki].faceLinks().set(facei)
                 && links_[bricki].faceLinks()[facei].T() != eye
                )
                {
                    aligned_ = false;
                }
            }
        }
    }

    // Set periodic face links

    setPeriodicFaceLinks();

    // Only set edge and vertex links if the bricks are structured. For
    // unstructured bricks, edge or vertex links may be ambiguous or undefined.

    if (structured_)
    {
        setEdgeAndVertexLinks();
    }
}

brickTopology::~brickTopology()
{}

labelList brickTopology::shortestFacePath
(
    const label& bricki,
    const label& brickj,
    const bool ignorePeriodicFaces
) const
{
    // Dijkstra's algorithm

    labelList skip(geo_.bricks().size(), 0);
    labelList dist(geo_.bricks().size(), 99999);
    labelList prev(geo_.bricks().size(), -1);

    dist[bricki] = 0;

    while (min(skip) == 0)
    {
        label d = 99999;
        label u = -1;

        forAll(skip, i)
        {
            if (!skip[i] && dist[i] < d)
            {
                u = i;
                d = dist[i];
            }
        }

        forAll(links_[u].faceLinks(), facei)
        if
        (
            links_[u].faceLinks().set(facei)
         && (
                !ignorePeriodicFaces
             || !links_[u].faceLinks()[facei].periodic()
            )
        )
        {
            const brickFaceLink& link = links_[u].faceLinks()[facei];

            label v = link.f1().parentBrick().num();

            if (!skip[v])
            {
                label d = dist[u] + 1;

                if (d < dist[v])
                {
                    dist[v] = d;
                    prev[v] = u;
                }
            }
        }

        skip[u] = 1;
    }

    labelList P;

    label u = brickj;

    if (prev[u] > -1)
    while (u > -1)
    {
        P.append(u);
        u = prev[u];
    }

    return reverseList(P);
}

}

}
