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

Briscola is developed by [Edo Frederix](mailto:edo.frederix@nrgpallas.com),
[Victor Habiyaremye](mailto:victor.habiyaremye@nrgpallas.com) and Gonzalo
Bonilla at NRG PALLAS, the Netherlands.

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
* OpenMPI (version 3 or higher)
* CMake (version 3.22 or higher)
* FFTW
* Eigen (optional)
* PETSc (optional)

From OpenFOAM, only the OpenFOAM library (libOpenFOAM.so) is linked. In turn,
this library links against Pstream (libPstream.so). So only the OpenFOAM and
Pstream libraries need to be compiled and discoverable from your environment.

OpenMPI and FFTW are required while Eigen and PETSc are optional. OpenMPI should
already be available through your OpenFOAM installation and is automatically
used from that. The FFTW, Eigen and PETSc package locations can be specified by
the `FFTW_HOME`, `EIGEN_HOME` and `PETSC_HOME` environment variables. If those
variables are not set, an attempt is done to find the respective packages in
system locations using the pkg-config tool. If that fails too, the compilation
process will complain for required packages that a they were not found. When
PETSc and/or Eigen are found, the linear system solvers offered by these
packages are compiled into Briscola. In turn, these can then be used as 1)
coarse grid solver in the multigrid solver or 2) as main solver using the Krylov
solver. Briscola's Krylov solver class is a wrapper to the solvers offered by
the PETSc framework.

If the FFTW or Eigen packages are not already on your system, they can be
installed with:

```
cd dependecies
./makeFFTW
./makeEigen
```

By default, these make scripts configure, compile and install their respective
packages to `$HOME/opt`. You can edit the make scripts if you want to specify
another location. After completion, both scripts will instruct you on which
environment variables to set. If you would like to use linear solvers from
PETSc, the PETSc package must be installed by yourself.

## Building Briscola

Unlike OpenFOAM, Briscola uses CMake for configuration, building and
installation. It is recommended to build Briscola using the CMake configure and
build *presets* called 'default' as defined in `CMakePresets.json`, i.e.:

```
cmake --preset default
cmake --build --preset default
```

This configures and builds the code with compiler optimization. Build files are
written to the `build` directory, which is automatically created. The build step
also automatically performs an installation, which by default is to your
`$FOAM_USER_LIBBIN` for libraries and `$FOAM_USER_APPBIN` for executables. In
this way, OpenFOAM's `wmake` behavior is fully mimicked. The presets
automatically use all parallel resources available on your machine.

It also possible to compile with debug flags. This is done using the 'debug'
presets replacing the 'default' ones. It is also possible to compile just a
single solver, by selecting specific targets. For example, to only compile the
`briscolaColocated` solver, do:

```
cmake --preset default
cmake --build --preset default --target briscolaColocated
```

Other custom targets include `libraries`, `solvers` and `utilities`. Note that
using specific targets does not automatically perform an installation anymore;
the built libraries and/or binaries will be in build/lib and build/bin.

To verify your Briscola installation, we provide a large number of unit test
applications. These can be built with the 'tests' preset (or 'tests-debug' to
disable compiler optimization). To build and run the tests, do:

```
cmake --preset default
cmake --build --preset tests
ctest --preset default
```

The `ctest` command reports the success rate of all performed tests. Note that
in default compilation mode (i.e., with compiler optimization) the build time of
some tests can take minutes.

Users are also free not to use CMake presets, if a more traditional use of CMake
is preferred. The CMake configuration accepts standard parameters like
`-DCMAKE_BUILD_TYPE=...` and `-DCMAKE_INSTALL_PREFIX=...` if needed.

## VSCode integration

CMake enables build and debug integration in editors like VSCode. For example,
using the 'CMake Tools' extension in VSCode the CMake presets can be used to
build from within VSCode. An important note is that since Briscola's CMake
configuration relies heavily on the OpenFOAM environment, it is key that VSCode
runs within the appropriate environment. When using SSH remotes, it is important
that OpenFOAM is loaded into your environment through `.bashrc` or `.cshrc`
files or similar. When multiple environments are needed on a remote machine, it
is recommended to use `code-server` instead, because SSH remotes handle
environments poorly (see this [this issue on
Github](https://github.com/microsoft/vscode-remote-release/issues/141)).

Briscola also supporst clangd integration, using the clangd VSCode extension. A
default `.clangd` configuration file is provided, which suppresses diagnostics
and so-called inlay hints. The reason for suppressing diagnostics is that most
source files will fail to compile because of a lack of appropriate environment.
As such, the clangd extension is only useful for definition tracking and inline
code suggestions.

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

## Development

Briscola is currently developed on Github. A simplified Gitflow model is used,
with the following rules:

* All developments are collected in the master branch.
* Changes must be made via a separate branch, named feature/\<description\>,
  bugfix/\<description\>, where \<description\> describes what was changed,
  e.g., feature/someNewScheme or bugfix/someSolver, etc.
* Once completed, a pull request must be created from this branch, and at least
  one reviewer must agree with the proposed changes before it can be merged.
* Users are encouraged to create forks and to propose changes via a pull request
* When adding new features, it is required to add test applications that verify
  the new feature implementation, and to add the feature to one of the cases or
  to create a new case that uses the feature.
* Once in a while a new release tag is created on the master branch. Briscola
  releases are in the format v\<major\>.\<minor\>.\<patch\>. Currently,
  development is mostly demand-driven, so versioning and releasing is ad hoc
  with little planning.

## Contact & Support

For bug reports or support, feel free to contact [Edo
Frederix](mailto:edo.frederix@nrgpallas.com).

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
