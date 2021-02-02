#!/bin/sh
set -ex
# Check if NTL is already installed
if [ ! -d "$HOME/ntl/lib" ] || [ -f "ntl-*.tar.gz" ] ; then
  wget https://libntl.org/ntl-11.4.3.tar.gz
  tar -xzvf ntl-*.tar.gz
  cd ntl-9.6.2/src && make && make install
else
  echo 'NTL is already installed';
  fi