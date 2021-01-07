# MSMPathfinder
Finding pathways from an Markov State Model.

This software package provides tools for finding pathways based on Markov state
models. It is required to represent the matrix as an row-normalized matrix,
where each row sums up to 1 and all entries are strictly greater 1.
This package contains a proof-of-concept implementation of an suggested
algorithm (see citations). If the matrix contains more than 100 highly
connected states or a few hundred sparsly connected states the MCMC is
preferable.

The essential functions are:
  - `MSMPathfinder paths`: see citation for further information
  - `MSMPathfinder mcmcsingle`: Markov chain Monte Carlo propagation

 For further information simply call both submodule with added `-h` flag.


# Documentation
All options are well documented and may be viewed by 'MSMPathfinder -h'.


# Citations
The underlying methods are based on the following article:
  - D. Nagel, A. Weber and G. Stock,
    *MSMPathfinder: Identification of pathways in Markov state models*,
    J. Chem. Theory Comput., 16, 7874 (2020);
    DOI: [10.1021/acs.jctc.0c00774](https://doi.org/10.1021/acs.jctc.0c00774)

We kindly ask you to cite this article if you use this software package for
published works.


# Licensing
This project was created by [moldyn-nagel](https://github.com/moldyn-nagel).

Copyright (c) 2020-2021, Daniel Nagel
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


# Installation
## Requirements
 required:
  -  **BOOST >= 1.49**
  -  **cmake >= 2.8**
  -  a **recent** C++ compiler (e.g. GNU g++ >= 4.9, must
     support C++11 standard)


## Quick-Start

To quickly get a working binary

```bash
# clone the repo
git clone https://github.com/moldyn/MSMPathfinder.git

# create a build folder inside the code directory
cd MSMPathfinder
mkdir build

# change to the build directory and prepare make
# if installition is not desired you can execute "cmake .." instead
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/my/installation/path

# then compile and install the package to /my/installation/path
# (or any other path you chose above) by invoking
make
make install
```
