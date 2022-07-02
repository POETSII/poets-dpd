#!/bin/bash

sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt install g++-11 g++ build-essential mpi-default-dev -y

if [[ ! -d lammps ]] ; then
    git clone https://github.com/lammps/lammps.git
fi
(
    cd lammps
    mkdir build_cpu
    (
        if [[ ! -x build_cpu/lmp ]] ; then
            cd build_cpu
            >&2 echo "Building CPU"
            /snap/bin/cmake -D PKG_DPD-BASIC=on -D PKG_OPENMP=on -D PKG_INTEL=on \
                -D PKG_OPT=on -D PKG_KOKKOS=on -DKokkos_ENABLE_OPENMP=yes \
                -D CMAKE_CXX_COMPILER=g++-11 -D CMAKE_CXX_FLAGS="-O3 -mavx2 -ffast-math" \
                ../cmake
            make -j$(nproc)
        fi
    )
    if [[ "$(ls /dev/nv*)" != "" ]] ; then
        if [[ ! -x build_gpu/lmp ]] ; then
            mkdir build_gpu
            (
                cd build_gpu
                cmake -D PKG_GPU=on -D PKG_DPD-BASIC=on -D PKG_OPENMP=on -D PKG_INTEL=on \
                    -D PKG_OPT=on -D PKG_KOKKOS=on  \
                    -D Kokkos_ENABLE_CUDA=yes -D Kokkos_ENABLE_SERIAL=yes  \
                    -D CMAKE_CXX_FLAGS="-O3 -mavx2 -ffast-math" ../cmake
                make -j$(nproc)
            )
        fi
        if [[ ! -x build_gpu_single/lmp ]] ; then
            mkdir build_gpu_single
            (
                cd build_gpu_single
                cmake -D PKG_GPU=on -D PKG_DPD-BASIC=on -D PKG_OPENMP=on \
                    -D PKG_OPT=on -D GPU_PREC=single \
                    -D CMAKE_CXX_FLAGS="-O3 -mavx2 -ffast-math" ../cmake
                make -j$(nproc)
            )
        fi

    fi
)
