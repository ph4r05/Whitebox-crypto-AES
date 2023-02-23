#!/bin/sh
set -ex
wget http://www.shoup.net/ntl/ntl-11.5.1
tar -xzvf ntl-11.5.1
cd ntl-11.5.1/src && ./configure && make && make install
