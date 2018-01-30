ECHMETCoreLibs
===

Introduction
---

ECHMETCoreLibs is a set of libraries that implements various computation methods to calculate chemical concentration equilibria in systems of general electrolytes. The result of the implemented calculations is primarily intended to serve as the base input for mathematical models of electromigration in capillary electrophoresis experimental setup.

ECHMETCoreLibs include a powerful equilibria solver that can solve concentration equilibria in systems containing arbitrary amount of general ampholytes. Additionaly, the solver can account for non-trivial complex-forming equilibira between arbitrary constituents in a system. The complex-forming equilibria model is based on the needs of affinity capillary electrophoresis (ACE) simulations and has limitations as to which components of a chemical system can form complexes. The limitations of the model are discussed further. (TODO)

The ECHMETCoreLibs project is divided into four separate libraries.

### libECHMETShared
Library with base components that are used by the rest of the libraries. This library also provides all physical constants and basic conversion functions.

### libSysComp
Task of `libSysComp` is to convert a chemical system description form the user to internal representation that is used for subsequent calculations. `libSysComp` also does all sanity checks of the input data.

### libCAES
`libCAES` implements the actual **C**omplex-forming and **A**cidobazic **E**quilibria **S**olver. It also provides conveniece functions to calculate numerical derivatives of individual ionic concentrations.

### libIonProps
This library implements calculations of certain ionic properties and non-ideality corrections. These corrections include Onsager-Fuoss correction of ionic mobilities and rudimentary correction for viscosity effects.

It shall be noted that these libraries are intended to be built and used together. Building the libraries individually or mixing libraries from different builds is untested and highly discouraged.

Building
---

Prior to building ECHMETCoreLibs make sure that you have all the necessary dependencies installed. ECHMETCoreLibs depend on the following tools and libraries:

 - [GNU Multiple Precision library (GMP)](https://gmplib.org/)
 - [GNU MPFR library](http://www.mpfr.org/)
 - [Eigen](http://eigen.tuxfamily.org/)
 - [CMake](https://cmake.org/)


### Linux/UNIX

In order to configure the project, proceed as follows

`cd` into the project's source directory
create a build directory (`mkdir build`) and `cd` into it
run `cmake -DCMAKE_BUILD_TYPE=Release`
run `make && make install`

If you do not have a system-wide installation of the Eigen library, you may specify a custom path with additional `-DEIGEN_INCLUDE_DIR=<path>` CMake parameter.

Similarly, if you do not have a system-wide installation of the GMP and MPFR libraries, you may specify custom paths with the following sequence of parameters: `-DMANUAL_HIPREC_LIBS_PATH=ON -DGMP_LIBRARY_BIN=<path_to_libgmp.so> -DMPFR_LIBRARY_BIN=<path_to_libmpfr.so> -DGMP_INCLUDE_DIR=<path_to_gmp.h> -DMPFR_INCLUDE_DIR=<path_to_mpfr.h>`

### Windows

To be added, you are on your own for now...

Examples
---

There is a sample reference tool included in the `testing` directory that demonstrates the intended use of the libraries and their API. The tool takes input from a JSON file. JSON files with some sample chemical systems can be found in `testing\testdata` directory. You may need to update the shell scripts in the `testing` directory to get the reference tool to build and run properly.

Licensing
---

The ECHMETCoreLibs project is distributed under the terms of **The GNU General Public License v3** (GNU GPLv3). See the enclosed `LICENSE` file for details.

As permitted by section 7. *Additional Terms* of The GNU GPLv3 license, the authors require that any derivative work based on ECHMETCoreLibs clearly refers to the origin of the software and its authors. Such reference must include the address of this source code repository (https://github.com/echmet/ECHMETCoreLibs) and names of all authors and their affiliation stated in section [Authors](#Authors) of this README file.

<a name="Authors"></a>
Authors
---

Michal Malý  
Pavel Dubský

Group of Electromigration and Chromatographic Methods (http://echmet.natur.cuni.cz)

Department of Physical and Macromolecular Chemistry  
Faculty of Science, Charles University, Czech Republic
