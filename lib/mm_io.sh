#! /bin/bash
#
cp mm_io.h /$HOME/include
#
g++-10 -c -Wall mm_io.c
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi