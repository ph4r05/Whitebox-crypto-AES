#!/bin/sh
set -e
# check to see if protobuf folder is empty
if [ ! -d "$HOME/ntl/lib" ]; then
  wget http://www.shoup.net/ntl/ntl-9.6.2.tar.gz
  tar -xzvf ntl-9.6.2.tar.gz
  cd ntl-9.6.2/src && ./configure PREFIX=$HOME/ntl && make && make install
else
  echo 'Using cached directory.';
fi