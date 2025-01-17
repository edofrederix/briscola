namespace Foam
{

namespace briscola
{

namespace fv
{

inline void truncatedTet::setVertexConstants()
{
    phi_.setSize(4);

    for (int i = 0; i < 4; i++)
        phi_[i] = (n_ & v_[i]) + C_;
}

inline truncatedTet::truncatedTet
(
    const vectorList& v,
    const vector& n,
    const scalar& C
)
:
    v_(v),
    n_(n),
    C_(C)
{
    setVertexConstants();
}

inline truncatedTet::truncatedTet(const truncatedTet& hx)
:
    v_(hx.v_),
    n_(hx.n_),
    C_(hx.C_),
    phi_(hx.phi_)
{}

inline scalar truncatedTet::volume() const
{
    label n = 0;
    forAll(phi_, i)
        n += phi_[i] >= 0;

    switch(n)
    {
        case 0:

            return 0;

        case 1:
        {
            vectorList p(3);
            vector q = Zero;
            label k = 0;

            for (int i = 0; i < 4; i++)
            if (phi_[i] >= 0)
            {
                q = v_[i];

                for (int j = 0; j < 4; j++)
                    if (i != j)
                        p[k++] = intersectionPoint(i,j);
            }

            return Foam::mag(tetVolume(p[0], p[1], p[2], q));
        }
        case 2:
        {
            vectorList p(4);
            labelList  q(4);
            labelList  s(2);

            label j = 0;
            label k = 0;

            for (int i = 0; i < 4; i++)
            if (phi_[i] >= 0)
            {
                s[j++] = i;

                for (int j = 0; j < 4; j++)
                if (i != j && phi_[j] < 0)
                {
                    p[k] = intersectionPoint(i,j);
                    q[k] = j;

                    k++;
                }
            }

            label a = (q[0] == q[2] ? 3 : 2);

            return
                Foam::mag(tetVolume(v_[s[0]], v_[s[1]], p[0], p[1]))
              + Foam::mag(tetVolume(v_[s[1]], p[2], p[3], p[0]))
              + Foam::mag(tetVolume(v_[s[1]], p[0], p[1], p[a]));
        }
        case 3:
        {
            vectorList p(3);
            vector q = Zero;
            label k = 0;

            for (int i = 0; i < 4; i++)
            if (phi_[i] < 0)
            {
                q = v_[i];

                for (int j = 0; j < 4; j++)
                    if (i != j)
                        p[k++] = intersectionPoint(i,j);
            }

            return
                Foam::mag(tetVolume(v_[0], v_[1], v_[2], v_[3]))
              - Foam::mag(tetVolume(p[0], p[1], p[2], q));
        }
        case 4:

            return Foam::mag(tetVolume(v_[0], v_[1], v_[2], v_[3]));

        default:

            FatalErrorInFunction
                << "Cannot compute truncated tet volume"
                << endl << abort(FatalError);

            return 0;
    }
}

}

}

}
