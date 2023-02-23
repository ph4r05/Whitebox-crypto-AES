#!/bin/sh
set -e
# check to see if protobuf folder is empty
# Lib `gmp` is required, configure include paths correctly, i.e.,
# ./configure PREFIX=$HOME/ntl GMP_PREFIX=$(brew --prefix)
if [ ! -d "$HOME/ntl/lib" ]; then
  wget https://libntl.org/ntl-11.5.1.tar.gz
  tar -xzvf ntl-11.5.1.tar.gz
  cd ntl-11.5.1/src && ./configure PREFIX=$HOME/ntl && make && make install
else
  echo 'Using cached directory.';
fi