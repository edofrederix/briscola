// Set coefficients and sources to 0 in IB
USys.A() *= (1.0 - IB.stagMask());
USys.b() *= (1.0 - IB.stagMask());

forAll(USys.A(), l)
{
    forAll(USys.A()[l], d)
    {
        forAllCells(USys.A()[l][d], i, j, k)
        {
            // Set central coefficients to 1 in IB
            USys.A()[l][d](i,j,k).center() += IB.stagMask()[l][d](i,j,k);

            // Modify stencils in IB-adjacent cells
            // We only apply the IBM on the finest mesh (l == 0)
            if ((l == 0) && (IB.stagWallAdjMask()[l][d](i,j,k) == 1))
            {
                if (IB.stagWallDist()[l][d](i,j,k).left() >= 0)
                {
                    scalar xi = IB.stagWallDist()[l][d](i,j,k).left();

                    scalar w1 = 2.0 - (2.0 - xi) * 8.0/3.0;
                    scalar w2 = -1.0 + (1.0 - xi) * 8.0/3.0;

                    if (xi < IB.xiStabilityFactor())
                    {
                        w1 = -2.0*xi/(1.0-xi);
                        w2 = xi/(2.0-xi);
                    }

                    scalar a0 = USys.A()[l][d](i,j,k).left();

                    USys.A()[l][d](i,j,k).left() = 0;

                    USys.A()[l][d](i,j,k).center() += a0*w1;

                    USys.A()[l][d](i,j,k).right() += a0*w2;
                }

                if (IB.stagWallDist()[l][d](i,j,k).bottom() >= 0)
                {
                    scalar xi = IB.stagWallDist()[l][d](i,j,k).bottom();

                    scalar w1 = 2.0 - (2.0 - xi) * 8.0/3.0;
                    scalar w2 = -1.0 + (1.0 - xi) * 8.0/3.0;

                    if (xi < IB.xiStabilityFactor())
                    {
                        w1 = -2.0*xi/(1.0-xi);
                        w2 = xi/(2.0-xi);
                    }

                    scalar a0 = USys.A()[l][d](i,j,k).bottom();

                    USys.A()[l][d](i,j,k).bottom() = 0;

                    USys.A()[l][d](i,j,k).center() += a0*w1;

                    USys.A()[l][d](i,j,k).top() += a0*w2;
                }

                if (IB.stagWallDist()[l][d](i,j,k).aft() >= 0)
                {
                    scalar xi = IB.stagWallDist()[l][d](i,j,k).aft();

                    scalar w1 = 2.0 - (2.0 - xi) * 8.0/3.0;
                    scalar w2 = -1.0 + (1.0 - xi) * 8.0/3.0;

                    if (xi < IB.xiStabilityFactor())
                    {
                        w1 = -2.0*xi/(1.0-xi);
                        w2 = xi/(2.0-xi);
                    }

                    scalar a0 = USys.A()[l][d](i,j,k).aft();

                    USys.A()[l][d](i,j,k).aft() = 0;

                    USys.A()[l][d](i,j,k).center() += a0*w1;

                    USys.A()[l][d](i,j,k).fore() += a0*w2;
                }

                if (IB.stagWallDist()[l][d](i,j,k).right() >= 0)
                {
                    scalar xi = IB.stagWallDist()[l][d](i,j,k).right();

                    scalar w1 = 2.0 - (2.0 - xi) * 8.0/3.0;
                    scalar w2 = -1.0 + (1.0 - xi) * 8.0/3.0;

                    if (xi < IB.xiStabilityFactor())
                    {
                        w1 = -2.0*xi/(1.0-xi);
                        w2 = xi/(2.0-xi);
                    }

                    scalar a0 = USys.A()[l][d](i,j,k).right();

                    USys.A()[l][d](i,j,k).right() = 0;

                    USys.A()[l][d](i,j,k).center() += a0*w1;

                    USys.A()[l][d](i,j,k).left() += a0*w2;
                }

                if (IB.stagWallDist()[l][d](i,j,k).top() >= 0)
                {
                    scalar xi = IB.stagWallDist()[l][d](i,j,k).top();

                    scalar w1 = 2.0 - (2.0 - xi) * 8.0/3.0;
                    scalar w2 = -1.0 + (1.0 - xi) * 8.0/3.0;

                    if (xi < IB.xiStabilityFactor())
                    {
                        w1 = -2.0*xi/(1.0-xi);
                        w2 = xi/(2.0-xi);
                    }

                    scalar a0 = USys.A()[l][d](i,j,k).top();

                    USys.A()[l][d](i,j,k).top() = 0;

                    USys.A()[l][d](i,j,k).center() += a0*w1;

                    USys.A()[l][d](i,j,k).bottom() += a0*w2;
                }

                if (IB.stagWallDist()[l][d](i,j,k).fore() >= 0)
                {
                    scalar xi = IB.stagWallDist()[l][d](i,j,k).fore();

                    scalar w1 = 2.0 - (2.0 - xi) * 8.0/3.0;
                    scalar w2 = -1.0 + (1.0 - xi) * 8.0/3.0;

                    if (xi < IB.xiStabilityFactor())
                    {
                        w1 = -2.0*xi/(1.0-xi);
                        w2 = xi/(2.0-xi);
                    }

                    scalar a0 = USys.A()[l][d](i,j,k).fore();

                    USys.A()[l][d](i,j,k).fore() = 0;

                    USys.A()[l][d](i,j,k).center() += a0*w1;

                    USys.A()[l][d](i,j,k).aft() += a0*w2;
                }
            }
        }
    }
}