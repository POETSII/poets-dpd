#!/bin/bash

sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt install g++-11 cmake g++ build-essential mpi-default-dev -y

if [[ ! -d lammps ]] ; then
    git clone https://github.com/lammps/lammps.git
fi
(
    cd lammps
    mkdir build_cpu
    (
        cd build_cpu
	>&2 echo "Building CPU"
        /snap/bin/cmake -D PKG_DPD-BASIC=on -D PKG_OPENMP=on -D PKG_INTEL=on \
            -D PKG_OPT=on -D PKG_KOKKOS=on -DKokkos_ENABLE_OPENMP=yes \
            -D CMAKE_CXX_COMPILER=g++-11 -D CMAKE_CXX_FLAGS="-O3 -mavx2 -ffast-math" \
            ../cmake
        make -j16
    )
    if [[ "$(ls /dev/nvidia*)" != "" ]] ; then
        mkdir build_gpu
        (
            cd build_gpu
		>&2 echo "Building GPU"
            /snap/bin/cmake -D PKG_GPU=on -D PKG_DPD-BASIC=on -D PKG_OPENMP=on -D PKG_INTEL=on \
                -D PKG_OPT=on -D PKG_KOKKOS=on  \
                -D Kokkos_ENABLE_CUDA=yes -D Kokkos_ENABLE_SERIAL=yes  \
                -D CMAKE_CXX_FLAGS="-O3 -mavx2 -ffast-math" ../cmake
            make -j4
        )
    fi
)
