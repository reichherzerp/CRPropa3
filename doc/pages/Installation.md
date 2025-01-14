# Installation

## Prerequisites
1. Install relevant packages
     ```sh
     brew update
     brew install python virtualenv hdf5 fftw cfitsio muparser libomp numpy swig llvm zlib
     ```

## Virtual Environment Setup
1. Create and activate the virtual environment:
     ```sh
     export VENV=crpropa
     python3 -m venv $(pwd)/$VENV
     source $(pwd)/$VENV/bin/activate
     pip install numpy
     ```
2. Verify Python in the virtual environment:
     ```sh
     which python
     ```
Output should be: $(pwd)/$VENV/bin/python
## Download CRPropa
1. Clone the CRPropa repository into the virtual environment directory:
     ```sh
     cd $(pwd)/$VENV
     git clone https://github.com/CRPropa/CRPropa3.git
     cd CRPropa3
     ```
## Build and Install CRPropa
1. Create a build directory:
    ```sh
    mkdir build && cd build
    ```
    
2. Run CMake:
    ```sh
    cmake .. \
      -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV \
      -DPYTHON_EXECUTABLE=$(which python) \
      -DFFTW_ROOT=/opt/homebrew \
      -DGSL_ROOT_DIR=/opt/homebrew \
      -DHEALPIX_ROOT=/opt/homebrew \
      -DENABLE_GALACTICMAGNETICLENS=FALSE
    ```
3. Compile and install:
    ```sh
    make -j$(sysctl -n hw.ncpu)
    make install
    make test
    ```


## Test installation
1. Test the installation:
    ```sh
    python -c "import crpropa; print(crpropa.pc)"
    ```
    
