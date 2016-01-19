#!/bin/sh

autoreconf --install

args="--prefix=/usr/local \
CPPFLAGS='-I/usr/local/include/flint' \
CFLAGS='-march=native -O3 -ffast-math -g -Wall -Wextra'"

echo
echo "----------------------------------------------------------------"
echo "Initialized build system. For a common configuration please run:"
echo "----------------------------------------------------------------"
echo
echo "./configure $args"
echo
