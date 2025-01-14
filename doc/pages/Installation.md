# Installation

## Prerequisites
1. Install relevant packages
     ```sh
     brew update
     brew install cmake gcc fftw gsl muparser hdf5
     pip install numpy
     ```

## Virtual Environment Setup
1. Create and activate the virtual environment:
     ```sh
     export VENV=crpropa
     python3 -m venv $(pwd)/venv/$VENV
     source $(pwd)/venv/$VENV/bin/activate
     ```
3. Verify Python in the virtual environment:
     ```sh
     which python
     # Output should be: $(pwd)/venv/$VENV/bin/python
     ```
## Download CRPropa
Clone the CRPropa repository into the virtual environment directory:
     ```sh
     cd $(pwd)/venv/$VENV
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
    crpropa_32 % python -c "import crpropa; print(crpropa.pc)"
    ```
    
