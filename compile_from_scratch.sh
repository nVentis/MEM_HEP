#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

compile_pkg ()
{
    cd $1
    mkdir -p build
    cd build
    cmake -DCMAKE_CXX_STANDARD=17 ..
    make install || cd ../.. && return 1
    cd ../..
}

cd source

local module_to_compile=""
for module_to_compile in $(ls .)
do
    compile_pkg $module_to_compile && echo "${GREEN}+++ Successfully compiled $module_to_compile +++${NC}" || echo "${RED}!!! Error [$?] while trying to compile $module_to_compile !!!${NC}" && cd .. && return 1
done

unset module_to_compile
unset RED
unset NC

cd ..
