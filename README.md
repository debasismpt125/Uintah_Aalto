# Uintah_Geotech
Uintah-MPM repository for geotechnical applications

# INSTALLATION

**1. Install pre-installation softwares**

sudo apt-get install subversion libhypre-dev petsc-dev \
libxml2-dev zlib1g-dev liblapack-dev cmake libglew-dev \
libxmu-dev g++ gfortran libboost-all-dev git \
libxrender-dev libxi-dev

**2. openmpi installation**

Follow the instructions on: https://sites.google.com/site/rangsiman1993/comp-env/program-install/install-openmpi

Example: ../configure --enable-mpi-thread-multiple --prefix=/usr/local

make

make install

**3. HYPRE installation**

Follow the instructions on: https://github.com/hypre-space/hypre

**4. visit Installation**

Follow the instructions on:  https://visit-dav.github.io/visit-website/releases-as-tables/

Example: sudo ./visit-install3_2_1 3.2.1 linux-x86_64-ubuntu20  /usr/local/visit

Add visit to path

export PATH="/usr/local/visit/bin:$PATH"

./visit

Or cd /usr/local/visit/bin

./visit

**5. PETSC Installation**

Follow the instructions on: https://petsc.org/release/install/download/

LD_LIBRARY_PATH=/usr/local/lib \
./configure --with-shared-libraries \
--with-debugging=O \
--with-mpi-dir=/usr/local\
--prefix=/home/jas/petsc

../src/configure --enable-debug --enable-all-components \
--with-boost=/usr --enable-wasatch_3p

**6. Compile  Uintah in computer**

Create directory named opt
cd to opt
../src/configure '--enable-optimize=-O3 -mfpmath=sse' --enable-mpm --without-fortran --with-mpi-lib=/usr/lib/x86_64-linux-gnu/openmpi/lib --with-mpi-include=/usr/lib/x86_64-linux-gnu/openmpi/include F77=gfortran
make

**Hypre installation (Optional)**

Follow Instructions on: https://github.com/hypre-space/hypre

cd hypre-2.18.2/src
./configure \
    --prefix=/usr/installed/hypre-2.18.2/gcc10.2.1-mpich3.4 \
    --enable-shared \
    --with-MPI-include=/usr/include/mpich \
    --with-MPI-lib-dirs=/usr/lib/mpich/lib \
    --with-MPI-libs='mpich' \
    CC=mpicc \
    CXX=mpicx

**Complile Uintah**

../src/configure '--enable-optimize=-O3 -mfpmath=sse' --enable-mpm  --enable-ice --without-fortran --with-mpi-lib=/usr/lib/x86_64-linux-gnu/openmpi/lib --with-mpi-include=/usr/lib/x86_64-linux-gnu/openmpi/include F77=gfortran  --with-hypre-lib=/home/debasis/Downloads/hypre-2.18.0/src/hypre/lib  --with-hypre-include=/home/debasis/Downloads/hypre-2.18.0/src/hypre/include





