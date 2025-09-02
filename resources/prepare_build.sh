#! /bin/bash

################################################################################
#
# This script configures the cmake environment
#
################################################################################

# get current directory
_cwd="$PWD"

# set dynadjust clone dir
_clone_dir="$_cwd"
# set dynadjust root dir
_root_dir="$_clone_dir/dynadjust"
# set build dir
_build_dir="$_root_dir/build_gcc"

echo " "
echo "Current directory: $_clone_dir"
echo "Root directory:    $_root_dir"
echo "Build directory:   $_build_dir"
echo " "

# 1. create build dir:
if [[ -d "$_build_dir" ]]; then
    echo "Cleaning out directory $_build_dir"
    cd "$_build_dir"
    rm -rf CMakeCache.txt CMakeFiles cmake_install.cmake dynadjust Makefile
    cd "$_cwd"
else
    echo "Creating new directory $_build_dir"
    mkdir "$_build_dir"
fi

cd "$_build_dir"

# 2. copy files:
echo "Copying Find...cmake files to build directory..."
#cp ../FindXercesC.cmake ./
#cp ../FindMKL.cmake ./
cp ../FindXSD.cmake ./

echo "Done."
echo " "

