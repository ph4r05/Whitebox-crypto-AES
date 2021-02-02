#!/bin/sh
set -e
# check to see if protobuf folder is empty
if [ ! -d "$HOME/ntl/lib" ] || [ -f "ntl-*.tar.gz" ] ; then
  wget https://libntl.org/ntl-11.4.3.tar.gz
  tar -xzvf ntl-*.tar.gz
  cd ntl-*/src && ./configure PREFIX="$HOME"/ntl && make && make install
else
  echo 'Using cached directory.';
fi