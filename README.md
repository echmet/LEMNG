LEMNG
===

Introduction
---
LEMNG (**L**inear **E**lectromigration **M**odel **N**ext **G**eneration) is a computer implementation of a mathematical model of linear theory of electromigration published at (TODO: Add articles when they get published). The model includes both acidobazic and complexation equilibria and is suitable for simulations of capillary zone electrophoresis experimental setups. LEMNG is built on top of [ECHMETCoreLibs](https://github.com/echmet/ECHMETCoreLibs) and uses its facilities to calculate equilibrium composition of solutions. LEMNG consists of a single library that can calculate mobilities, composition and shapes of all eigenzones for a system of given analytical composition and physical-chemical parameters. It can plot the expected electrophoregrams using HVL-R function profiles as well. In terms of complex-forming equilibria, LEMNG is constrained by the same limitations as ECHMETCoreLibs. Please refer to ECHMETCoreLibs documentation for further details.

Building
---
LEMNG requires the following tools to be built:

- C++14-aware compiler
- [Eigen matrix library](https://eigen.tuxfamily.org/)
- [CMake build system](http://cmake.org)
- [ECHMETCoreLibs](https://github.com/echmet/ECHMETCoreLibs)

### Generate makefiles with CMake
- Linux
1. `cd` into the source directory
2. Run `mkdir build` and `cd build`
3. Run `cmake .. -DCMAKE_BUILD_TYPE=Release -DEIGEN_INCLUDE_DIR=<path_to_eigen_library> -DECHMET_CORE_LIBS_DIR=<path_to_ECHMETCoreLibs_installation>`
4. Run `make` and `make install`
- Windows
1. You may use the CMake GUI tool the set up the required configuration options as it is described in the Linux section
2. Generate the appropriate project files. MSVC2015 is currently the recommended compiler on Windows although MinGW may be used as well.

Examples
---
There is a reference implementation available in `ref_tool/ref_tool.cpp`. To simplify the building process there is a series of shell scripts available. In order to build the reference tool, set the required paths accordingly to your setup in `ref_tool_glob.sh` and run `build_ref_tool.sh`. Project for MSVC is not available at the moment.

Licensing
---
The LEMNG project is distributed under the terms of **The GNU General Public License v3** (GNU GPLv3). See the enclosed `LICENSE` file for details.

As permitted by section 7. *Additional Terms* of The GNU GPLv3 license, the authors require that any derivative work based on LEMNG clearly refers to the origin of the software and its authors. Such reference must include the address of this source code repository (https://github.com/echmet/LEMNG) and names of all authors and their affiliation stated in section [Authors](#Authors) of this README file.

Referenceable releases
---
- LEMNG 0.6.0 [![DOI](https://zenodo.org/badge/113085202.svg)](https://zenodo.org/badge/latestdoi/113085202)

<a name="Authors"></a>
Authors
---
Michal Malý  
Pavel Dubský

Group of Electromigration and Chromatographic Methods (http://echmet.natur.cuni.cz)

Department of Physical and Macromolecular Chemistry  
Faculty of Science, Charles University, Czech Republic

