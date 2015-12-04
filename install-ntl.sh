#!/bin/sh
set -ex
wget http://www.shoup.net/ntl/ntl-9.6.2.tar.gz
tar -xzvf ntl-9.6.2.tar.gz
cd ntl-9.6.2/src && make && make install
