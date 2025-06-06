# About Briscola

## About

Briscola is a BRIck-structured Staggered and COLlocAted CFD code, aimed at
parallel high performance and high fidelity PDE simulation in relatively simple
but real-world geometries. The goal of Briscola is to be lean, efficient but
generic, to allow for high fidelity (e.g., LES or DNS) simulation of realistic
problems, within reasonable parallel computational effort.

Key to Briscola is a 'brick-structured' approach, which allows for the use of
fast and efficient geometric multigrid solvers. Parallel domain decomposition is
designed in such a way that the geometric multigrid solvers remain efficient
while keeping parallel communication at a minimum. Briscola uses OpenFOAM for
primitive data types, simple IO, parallel communication and compilation. Note
that from OpenFOAM, Briscola requires only the libOpenFOAM, libPstream and
libOSspecific header files and libraries.

## Developers

Briscola is developed by [Edo Frederix](mailto:edofrederix@gmail.com),
[Victor Habiyaremye](mailto:habiyaremye@nrg.eu) and
[Gonzalo Bonilla](mailto:bonillamoreno@nrg.eu) at NRG, the Netherlands.

## License

Briscola is published under the GNU GPL Version 3 license.

Briscola is distributed under the European Dual Use Codification N: EU DuC=N.
Goods labeled with an EU DuC (European Dual-Use Codification) not equal to 'N'
are subject to European and national export authorization when exported from the
EU and may be subject to national export authorization when exported to another
EU country as well. Even without an EU DuC, or with EU DuC 'N', authorization
may be required due to the final destination and purpose for which the goods are
to be used. No rights may be derived from the specified EU DuC or absence of an
EU DuC.

## Dependencies

Briscola depends on the following third-party packages:

* OpenFOAM (foundation version 12)
* OpenMPI (version 3 or later)
* FFTW
* Eigen
* PETSc
* SuperLU (optional)

From OpenFOAM, only the OpenFOAM library (libOpenFOAM.so) is linked. In turn,
this library links against Pstream (libPstream.so). So only the OpenFOAM and
Pstream libraries need to be compiled and discoverable from your environment

OpenMPI, FFTW, Eigen and PETSc are required. OpenMPI should already be available
through your OpenFOAM installation and is automatically used from that. The
FFTW, Eigen and PETSc package locations can be specified by the `FFTW_HOME`,
`EIGEN_HOME` and `PETSC_HOME` environment variables. If those variables are not
set, an attempt is done to find the respective packages in system locations
using the pkg-config tool. If that fails too, the compilation process will
complain that a package was not found. If these packages are not already on your
system, FFTW and Eigen can be installed with:

```
cd dependecies
./makeFFTW
./makeEigen
```

By default, these make scripts configure, compile and install their respective
packages to `$HOME/opt`. You can edit the make scripts if you want to specify
another location. Both scripts will instruct you on which environment variables
to set. The PETSc package must be installed by yourself.

The sparse direct solvers of SuperLU and SuperLU_DIST are used if these packages
are available, but they are not required nor important. The compiler checks for
the existence of the `SUPERLU_HOME` and `SUPERLU_DIST_HOME` environment
variables. If they exist, Briscola's SuperLU and SuperLU_DIST-based solvers are
compiled. It is interfaced via the Eigen and PETSc support functions. The
SuperLU and SuperLU_DIST packages can optionally be built from the
`dependencies` directory too with the `makeSuperLU` and `makeSuperLU_DIST`
scripts.

## Building Briscola

Once the dependencies above are met, Briscola can be built with

```
./Allwmake
```

This script uses OpenFOAM's wmake build system. For faster compilation, make
sure that the `WM_NCOMPPROCS` environment variable is set equal to the number of
processors available on your machine for compilation.

## Documentation

The code is self-documented via Doxygen. The Doxygen output can be generated
with

```
doxygen doc/Doxyfile
```

This requires doxygen, perl and graphviz to be installed. HTML output is written
to doc/Doxygen/html and can be best viewed by opening
doc/Doxygen/html/index.html.

There is also a short high-level documentation available that can be viewed
[HERE](doc/chapters/0_start.md).

## Contact & Support

For bug reports or support, feel free to contact [Edo
Frederix](mailto:edofrederix@gmail.com).

## Disclaimer

Briscola is provided by the copyright holders and contributors "as-is" and any
express or implied warranties, including, but not limited to, the implied
warranties of merchantability and fitness for a particular purpose are
disclaimed. In no event shall the copyright owner or contributors be liable for
any direct, indirect, incidental, special, exemplary, or consequential damages
(including, but not limited to, procurement of substitute goods or services;
loss of use, data, or profits; or business interruption) however caused and on
any theory of liability, whether in contract, strict liability, or tort
(including negligence or otherwise) arising in any way out of the use of this
software, even if advised of the possibility of such damage.
