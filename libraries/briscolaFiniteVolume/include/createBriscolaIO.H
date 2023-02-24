IO io(fvMsh);

if (runTime.value() != 0.0)
{
    io.read<colocated>();

    if (fvMsh.topology().structured())
        io.read<staggered>();
}
