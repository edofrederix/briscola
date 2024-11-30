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

* OpenFOAM (only libOpenFOAM, libPstream and libOSspecific)
* OpenMPI (version 3 or later)
* FFTW3
* Eigen3
* Intel oneAPI MKL (optional)
* SuiteSparse (optional)
* SuperLU (optional)

FFTW3 and Eigen3 are required, and their locations should be specified by the
FFTW_HOME and EIGEN_HOME environment variables.

The sparse direct solvers of MKL, SuiteSparse and SuperLU are used if these
packages are available. The compiler checks for the existence of the MKLROOT,
SUITESPARSE_HOME and SUPERLU_HOME environment variables. If they exist, the
Pardiso, UmfPack and SuperLU solvers are compiled, respectively. They are
interfaced via the Eigen support functions.

## Documentation

The code is self-documented via Doxygen. The Doxygen output can be generated
with

```
doxygen doc/Doxyfile
```

This requires doxygen and graphviz to be installed. There is also a short
high-level documentation available that can be viewed
[HERE](doc/chapters/start.md).

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
