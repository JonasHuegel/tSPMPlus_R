# tSPM+ R package
This packages provides an R interface to the C++ tSPM+ library.
For the example use cases please read the vignettes in vignette folder.

## Hardware Requirements
tSPM+ requires a usual computer with a sufficient amount of memory (we recommend at least 8GB).

## Software Requirements
### OS Requirements
This R packages is tested  on Linux and Windows. It should also compile on macOS as long the user uses a compiler that enables openMP. We tested it on:
* Windows 10.0.19045
* Ubuntu 22.10

### Depencies to build the package:
* R 4.2.1+
* R-tools (windows only) version corresponding to R
* openMP
* C++ compiler supporting C++17 and openMP (we recommend gcc

#### R-packages depencies
* RcppParallel
* Rccp (>= 1.0.10)
* data.table
* dplyr
* utils
* memuse
* devtools

## Installation
This package wraps the tSPM+ C++ library, which is included as a git submodule. Adding the `--recursive` argumente when installing the package is not anymore required.
Install in R with:
`devtools::install_github("jonashuegel/tSPMPlus_R")`

## Usage
You can find two vignettes with example use cases in the vignette folder.

## Citation
Please cite our publication when using the package: Hügel, J., Sax, U., Murphy, S. N. & Estiri, H. tSPM+; a high-performance algorithm for mining transitive sequential patterns from clinical data. arXiv [cs.LG] (2023) doi:10.48550/arXiv.2309.05671

## License
This work is under the MIT License.
