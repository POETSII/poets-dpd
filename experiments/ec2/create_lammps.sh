#!/bin/bash

sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt install g++-11

git clone https://github.com/lammps/lammps.git
(
    cd lammps
    mkdir build_cpu
    (
        cd build_cpu
        cmake -D PKG_GPU=on -D PKG_DPD-BASIC=on -D PKG_OPENMP=on -D PKG_INTEL=on \
            -D PKG_OPT=on -D PKG_KOKKOS=on -DKokkos_ENABLE_OPENMP=yes \
            -D CMAKE_CXX_COMPILER=g++-11 -D CMAKE_CXX_FLAGS="-O3 -mavx2 -ffast-math" \
            ../cmake
    )
    if [[ "$(ls /dev/nv*)" != "" ]] ; then
        mdkir build_gpu
        (cd build_gpu
            cmake -D PKG_GPU=on -D PKG_DPD-BASIC=on -D PKG_OPENMP=on -D PKG_INTEL=on \
                -D PKG_OPT=on -D PKG_KOKKOS=on -DKokkos_ENABLE_OPENMP=no \
                -D Kokkos_ENABLE_CUDA=yes -D Kokkos_ENABLE_SERIAL=yes  \
                -D CMAKE_CXX_FLAGS="-O3 -mavx2 -ffast-math" ../cmake
        )
    fi
)
