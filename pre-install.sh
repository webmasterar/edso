#! /bin/sh

tar -xvf vcflib.tar.gz
cd vcflib/tabixpp/htslib
make
cd ../..
make
cd ..
