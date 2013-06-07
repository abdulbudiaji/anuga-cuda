#!/bin/bash
rm -Rf build
rm ./siteconf.py
./configure.py --cuda-root=$CUDA_ROOT --prefix=user --no-use-shipped-boost --boost-inc-dir=/apps/boost/1.46.1/include/ --boost-lib-dir=/apps/boost/1.46.1/lib/ --cuda-inc-dir=/apps/cuda/4.2.9/include/ --boost-compiler=gcc44
