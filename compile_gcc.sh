#!/bin/bash

OS=`uname`
echo $OS

CC=gcc-5
CXX=g++-5

FRAMEWORK=/usr/local/Cellar/r/3.2.0_1/ 

  
  
SRC=lpm3.cpp
SHAREDOBJECT=lpm3.so 
OBJECT=lpm3.o

rm -f $SHAREDOBJECT $OBJECT 
$CXX -I/usr/local/Cellar/r/3.2.0_1/R.framework/Resources/include -DNDEBUG  -I/usr/local/opt/gettext/include -I/usr/local/opt/readline/include  -I/usr/local/include  -fopenmp -fPIC  -g -O2  -c $SRC -o $OBJECT 
$CXX -dynamiclib -lgomp -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/usr/local/Cellar/r/3.2.0_1/R.framework/Resources/lib -L/usr/local/opt/gettext/lib -L/usr/local/opt/readline/lib -L/usr/local/lib -o $SHAREDOBJECT $OBJECT -F$FRAMEWORK -framework R -lintl -Wl,-framework -Wl,CoreFoundation
  
  
  
