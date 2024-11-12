#!/bin/bash


echo 'installing coconut and curl.h'

# check for GCC

echo 'checking for GCC'

if ! [ -x "$(command -v gcc)" ]; then
  echo 'error: gcc is not installed. install it an try again :)' >&2
  exit 1
fi

# install dependencies locally
mkdir -p $HOME/local
cd $HOME/local

# download and install libcurl
echo 'downloding curl'

wget https://curl.se/download/curl-7.88.1.tar.gz
tar -xzvf curl-7.88.1.tar.gz
cd curl-7.88.1
./configure --prefix=$HOME/local
make && make install

# compile
echo 'compiling'

gcc -o coconut /path/to/coconut.c -I$HOME/local/include -L$HOME/local/lib -lcurl

echo 'Installation complete.'
