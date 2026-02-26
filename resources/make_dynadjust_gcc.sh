#! /bin/bash

################################################################################
#
# This script clones the latest version, configures the environment, then
# builds and installs DynAdjust using:
#   gcc 10.1.1
#   boost 1.69.0
#   xerces-c 3.1.4
#
################################################################################

# set defaults
_script="make_dynadjust_gcc.sh"
_auto=0 # default option is to ask for user input
_debug=0 # default option is to build release variant
_clone=0 # default option is to clone afresh
_test=0 # default option is to skip cmake tests
_install=0 # default option is to install binaries
_binary="all" # default option is to build all binaries
# _binaries is an array of binary names used to select which binary to build based
# on user input. The "...wrapper" binaries are ordered first, so that if the user
# enters part of a name (e.g. "imp"), the script will select the first occurrence
# of the binary wrapper (e.g. "dnaimportwrapper") rather than the shared object.
_binaries=(dnaimportwrapper dnaimport dnageoidwrapper dnageoid dnareftranwrapper dnareftran dnasegmentwrapper dnasegment dnaadjustwrapper dnaadjust dnaplotwrapper dnaplot dynadjust)

#################################################################################
# Common functions
#
# display example message
function example {
    echo -e "examples:"
	echo -e "  $_script --auto --no-clone --test --no-install"
}

# display usage message
function usage {
    echo -e "usage: $_script [options]\n"
}

# display help message (calls usage and example)
function help {
    echo ""
    usage
    echo -e "options:"
    echo -e "  -a [ --auto ]        Run automatically with no user interaction."
    echo -e "  -b [ --binary ] arg  Build a specific binary (e.g. \"dnaimport\" or \"dnaadjustwrapper\")."
    echo -e "                       By default, \"all\" binaries are built."
    echo -e "  -c [ --no-clone ]    By default, the latest version will be cloned from GitHub"
	echo -e "                       into the current directory, using:"
	echo -e "                         git clone https://github.com/geoscienceaustralia/DynAdjust.git"
    echo -e "                       Provide this option if building source from a local copy, e.g.:"
	echo -e "                         $ wget https://github.com/geoscienceaustralia/DynAdjust/archive/refs/tags/v1.1.0.tar.gz -O DynAdjust-1.1.0.tar.gz"
    echo -e "                         $ tar xzvf DynAdjust-1.1.0.tar.gz"
    echo -e "                         $ cd DynAdjust-1.1.0/"
	echo -e "                         $ bash ./resources/make_dynadjust_gc.sh (this script)"
    echo -e "  -d [ --debug ]       Compile debug version."
    echo -e "  -n [ --no-install ]  Do not install binaries."
    echo -e "  -t [ --test ]        Run cmake tests."
    echo -e "  -h [ --help ]        Prints this help message.\n"
    example
    echo ""
}
#################################################################################

# get argument parameters
while [[ "$1" != "" ]];
do
   case $1 in
   -a  | --auto ) 	
   				  	_auto=1 # run automatically (or silently)
			    	;;
   -b  | --binary ) shift&&_binary=$1 # get the name of the binary to be built
			    	;;
   -d  | --debug )  
   					_debug=1 # compile debug variant
			        ;;
   -c  | --no-clone )
   					_clone=1 # do not clone from GitHub
			        ;;
   -t  | --test )   
   					_test=1 # run tests
			        ;;
   -n  | --no-install )
   					_install=1 # do not install binaries
			        ;;
   -h   | --help )  help
                    exit
                    ;;
   *)                     
                    echo -e "\n$_script: illegal option $1"
                    help
					exit 1 # error
                    ;;
    esac
	shift
done

# opt installation folder
OPT_DYNADJUST_PATH=/opt/dynadjust
OPT_DYNADJUST_GCC_PATH=/opt/dynadjust/gcc
DYNADJUST_INSTALL_PATH=/opt/dynadjust/gcc/1_3_0

# version info
_version="1.3.0"

echo -e "\n==========================================================================="
echo -e "DynAdjust $_version build configuration options..."

if [[ $_debug -eq 1 || $_test -eq 1 ]]; then
	echo -e " - debug variant."
else
	echo -e " - release variant."
fi

if [[ $_binary = "all" ]]; then
	echo -e " - building all binaries."
else
	# find equivalence of the whole word
	#if echo "${_binaries[@]}" | fgrep --word-regexp "$_binary"; then	
	
	# find equivalence with any part of a word
	# this allows dnaimportwrapper to be selected if imp is entered on 
	# the command line
	if [[ ${_binaries[@]} =~ ${_binary} ]]; then
		for binary in ${_binaries[@]}; do
			if [[ ${binary} =~ "$_binary" ]]; then
				_binary=${binary}
				echo -e " - building $_binary only."
				break
			fi
		done
	else
		echo -e " - building all binaries (I don't know what $_binary is)."
		_binary="all"
	fi
fi

if [[ $_auto -eq 1 ]]; then
	echo -e " - build automatically (no user input)."
else
	echo -e " - run interactively (ask for user input)."
fi

# get current directory
_cwd="$PWD"

if [[ $_clone -eq 1 ]]; then
	echo -e " - do not clone, but build from local source."
	# As per help, the user is expected to be in this directory once
	# a local copy has been extracted.
	_clone_dir="$_cwd"

else
	echo -e " - clone a fresh copy from GitHub."
	# set directory to which git will clone the latest
	# DynAdjust repo
	_clone_dir="$_cwd/DynAdjust"
fi

if [[ $_test -eq 1 ]]; then
	echo -e " - run tests."
fi

if [[ $_install -eq 1 ]]; then
	echo -e " - do not install."
else
	echo -e " - install binaries to ${DYNADJUST_INSTALL_PATH}."
fi

# set dynadjust root dir
_root_dir="$_clone_dir/dynadjust"
# set dynadjust root dir
_test_dir="$_clone_dir/sampleData"
# set build dir
_build_dir="$_root_dir/build_gcc"
# set clone url
_clone_url="https://github.com/geoscienceaustralia/DynAdjust.git"

# usr bin directory
BIN_FOLDER="~/bin"
eval BIN_FOLDER_FULLPATH="$BIN_FOLDER"

echo -e "\n==========================================================================="
echo -e "\nBuild and installation of DynAdjust $_version...\n"
if [[ $_clone -eq 0 ]]; then
	echo "Repository settings:"
	echo "  Git repo:      $_clone_url"
fi
echo "Build settings:"
echo "  Current dir:   $_cwd"
if [[ $_clone -eq 0 ]]; then
	echo "  Clone dir:     $_clone_dir"
fi
echo "  Build dir:     $_build_dir"

if [[ $_install -eq 0 ]]; then
	echo "Installation settings:"
	echo "  Install dir:   $DYNADJUST_INSTALL_PATH"
	echo "  User bin dir:  $BIN_FOLDER_FULLPATH"
fi

if [[ $_test -eq 1 ]]; then
	echo "Test settings:"
	echo "  Test dir:      $_test_dir"
fi

#
# determine whether user needs prompting
case ${_auto} in
    0) # perform interactive build
        echo " "
		read -r -p "Is this ok [Y/n]: " response;;
    *) # build without asking
        response="y";;
esac

if [[ "$response" =~ ^([nN][oO]|[nN])$ ]]
then    
    exit
else
    echo " "
fi

# INSTALL DYNADJUST
# 1. clone from GitHub (if required):
if [[ $_clone -eq 0 ]]; then
	echo " "
	echo "Cloning DynAdjust..."
	git clone "$_clone_url" || echo -e "Okay, let's assume we already have a previously cloned version.\n"
fi

# check if root directory exists
if [[ -d "$_root_dir" ]]; then
    # Good, directory exists
	echo -e "\nSuccessfully found root directory: $_root_dir\n"
else
    echo -e "\nerror: can't find root directory: $_root_dir\n"
	echo -e "  If building from a local copy, please ensure this script is executed"
	echo -e "  from the parent directory containing the following directories:"
	echo -e "    ./dynadjust"
	echo -e "    ./resources"
	echo -e "    ./sampleData\n"
	echo -e "  For example:"
	echo -e "    $ tar xzvf DynAdjust-1.1.0.tar.gz"
	echo -e "    $ cd DynAdjust-1.1.0/"
	echo -e "    $ ./resources/make_dynadjust_gcc.sh\n"
	exit 1 # error
fi

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

# 3. copy cmake files:
echo "Copying Find*.cmake files to build directory..."
#cp ../FindXercesC.cmake ./
#cp ../FindMKL.cmake ./
cp ../FindXSD.cmake ./

# check if MKLROOT environment variable has been set, which means
# intel MKL has been installed.
if [ -z "$MKLROOT" ]
then
    # If not, find setvars.h and execute to set environment variables
	if [ -e /opt/intel/oneapi/setvars.sh ]
	then
		echo "Setting environment variables for intel..."
		source /opt/intel/oneapi/setvars.sh
	fi
else
	echo "Intel root: $MKLROOT"
fi

REL_BUILD_TYPE="Release"
DBG_BUILD_TYPE="Debug"
THIS_BUILD_TYPE=$REL_BUILD_TYPE

# 4. build:
# Force build type to Debug for --debug or --test options
if [[ $_debug -eq 1 || $_test -eq 1 ]]; then
	# debug
	THIS_BUILD_TYPE=$DBG_BUILD_TYPE
fi

echo " "
echo "Building DynaNet ($THIS_BUILD_TYPE)..."
echo " "

gcc_version=$(gcc -v 2>&1 | tail -1 | awk '{print $1 " " $2 " " $3}')

echo "$gcc_version"
echo " "

#
# determine whether to prepare cmake files with testing or not
case ${_test} in
    0) # skip tests
        echo -e "cmake -DBUILD_TESTING=OFF -DUSE_MKL=ON -DCMAKE_BUILD_TYPE=$THIS_BUILD_TYPE ..\n"
		cmake -DBUILD_TESTING="OFF" -DUSE_MKL=ON -DCMAKE_BUILD_TYPE="$THIS_BUILD_TYPE" .. || exit 1;;
    *) # run cmake tests with code coverage
        echo -e "cmake -DBUILD_TESTING=ON -DUSE_MKL=ON -DCMAKE_BUILD_TYPE=$THIS_BUILD_TYPE ..\n"
		cmake -DBUILD_TESTING="ON" -DUSE_MKL=ON -DCMAKE_BUILD_TYPE="$THIS_BUILD_TYPE" .. || exit 1;;
esac

echo -e "\n==========================================================================="


# Make!
if [[ $_binary = "all" ]]; then
	echo -e "Building DynAdjust $_version...\n"
	make -j $(nproc) || exit 1
else
	echo -e "Building DynAdjust (${_binary}) $_version...\n"
	make -j $(nproc) ${_binary} || exit 1
fi

echo " "

case ${_test} in
    1) # run cmake tests
        echo -e "==========================================================================="
        echo -e "Testing DynAdjust $_version...\n"
        make CTEST_OUTPUT_ON_FAILURE=1 test
        # get current directory
        pwd
	# capture coverage info
	lcov --directory . --capture --output-file lcov.info
	# filter out system and test code
	lcov --remove lcov.info 'tests/*' '/usr/*' --output-file lcov.info
	# list before upload
	lcov --list lcov.info 
	;;
esac

#
# determine if user needs prompting
case ${_auto} in
    0)  # Interactive mode (check if the user has 
	    # already decided not to install)
		if [[ $_install -eq 1 ]]; then
			# user did not want to install, so don't ask
			optresponse="n"
		else
			echo " "
			read -r -p "Install DynAdjust to $DYNADJUST_INSTALL_PATH [Y/n]: " optresponse
		fi
        ;;
    *)  # proceed without asking. Assume installation
        optresponse="y"
		
		# user did not want to install, so alter the 
		# response accordingly
		if [[ $_install -eq 1 ]]; then
			optresponse="n"
		fi
		;;
esac

_lib_ext="so"

if [[ "$optresponse" =~ ^([nN][oO]|[nN])$ ]]
then    
    echo " "
else

	# create install dirs:
	if [[ $_install -eq 0 ]]; then
		if [[ ! -d "$BIN_FOLDER_FULLPATH" ]]; then
			echo " "
			echo "Making $BIN_FOLDER_FULLPATH"
			mkdir "$BIN_FOLDER_FULLPATH"
		fi
	fi

    # 1 GET OS, DISTRO AND TOOLS
	if [[ "$OSTYPE" == "darwin"* ]]; then
		# Mac OSX
		_lib_ext="dylib"
	fi

	if [[ ! -d $OPT_DYNADJUST_PATH ]]; then
		sudo mkdir $OPT_DYNADJUST_PATH
	fi

	if [[ ! -d $OPT_DYNADJUST_GCC_PATH ]]; then
		sudo mkdir $OPT_DYNADJUST_GCC_PATH
	fi

	if [[ ! -d $DYNADJUST_INSTALL_PATH ]]; then
		sudo mkdir $DYNADJUST_INSTALL_PATH
	fi

	echo "Copying libraries and binaries to $DYNADJUST_INSTALL_PATH ..."

	if [[ ( -e "../../bin/dynadjust" ) && ( "all" == "$_binary" || "$_binary" =~ "dynadjust" ) ]]; then
		sudo cp ../../bin/dynadjust "$DYNADJUST_INSTALL_PATH/"
		ln -sf "$DYNADJUST_INSTALL_PATH/dynadjust" "$BIN_FOLDER_FULLPATH/dynadjust"
		sudo ln -sf "$DYNADJUST_INSTALL_PATH/dynadjust" /opt/dynadjust/dynadjust
		echo " - dynadjust"
	fi

	if [[ ( -e "../../bin/dnaadjust" ) && ( "all" == "$_binary" || "$_binary" =~ "adjust" ) ]]; then
		sudo cp ../../bin/libdnaadjust.$_lib_ext "$DYNADJUST_INSTALL_PATH/"
		sudo cp ../../bin/dnaadjust "$DYNADJUST_INSTALL_PATH/"
		ln -sf "$DYNADJUST_INSTALL_PATH/dnaadjust" "$BIN_FOLDER_FULLPATH/dnaadjust"
		ln -sf "$DYNADJUST_INSTALL_PATH/libdnaadjust.$_lib_ext"  "$BIN_FOLDER_FULLPATH/libdnaadjust.$_lib_ext"
		sudo ln -sf "$DYNADJUST_INSTALL_PATH/libdnaadjust.$_lib_ext" /opt/dynadjust/libdnaadjust.$_lib_ext 
		sudo ln -sf "$DYNADJUST_INSTALL_PATH/dnaadjust" /opt/dynadjust/dnaadjust
		echo " - dnaadjust, libdnaadjust.$_lib_ext"
	fi

	if [[ ( -e "../../bin/dnaimport" ) && ( "all" == "$_binary" || "$_binary" =~ "import" ) ]]; then
		sudo cp ../../bin/libdnaimport.$_lib_ext "$DYNADJUST_INSTALL_PATH/"
		sudo cp ../../bin/dnaimport "$DYNADJUST_INSTALL_PATH/"
		ln -sf "$DYNADJUST_INSTALL_PATH/dnaimport" "$BIN_FOLDER_FULLPATH/dnaimport"
		ln -sf "$DYNADJUST_INSTALL_PATH/libdnaimport.$_lib_ext"  "$BIN_FOLDER_FULLPATH/libdnaimport.$_lib_ext"
		sudo ln -sf "$DYNADJUST_INSTALL_PATH/libdnaimport.$_lib_ext" /opt/dynadjust/libdnaimport.$_lib_ext
		sudo ln -sf "$DYNADJUST_INSTALL_PATH/dnaimport" /opt/dynadjust/dnaimport
		echo " - dnaimport, libdnaimport.$_lib_ext"
	fi

	if [[ ( -e "../../bin/dnareftran" ) && ( "all" == "$_binary" || "$_binary" =~ "reftran" ) ]]; then
		sudo cp ../../bin/libdnareftran.$_lib_ext "$DYNADJUST_INSTALL_PATH/"
		sudo cp ../../bin/dnareftran "$DYNADJUST_INSTALL_PATH/"
		ln -sf "$DYNADJUST_INSTALL_PATH/dnareftran" "$BIN_FOLDER_FULLPATH/dnareftran"
		ln -sf "$DYNADJUST_INSTALL_PATH/libdnareftran.$_lib_ext" "$BIN_FOLDER_FULLPATH/libdnareftran.$_lib_ext"
		sudo ln -sf "$DYNADJUST_INSTALL_PATH/libdnareftran.$_lib_ext" /opt/dynadjust/libdnareftran.$_lib_ext
		sudo ln -sf "$DYNADJUST_INSTALL_PATH/dnareftran" /opt/dynadjust/dnareftran
		echo " - dnareftran, libdnareftran.$_lib_ext"
	fi

	if [[ ( -e "../../bin/dnageoid" ) && ( "all" == "$_binary" || "$_binary" =~ "geoid" ) ]]; then
		sudo cp ../../bin/libdnageoid.$_lib_ext "$DYNADJUST_INSTALL_PATH/"
		sudo cp ../../bin/dnageoid "$DYNADJUST_INSTALL_PATH/"
		ln -sf "$DYNADJUST_INSTALL_PATH/dnageoid" "$BIN_FOLDER_FULLPATH/dnageoid"
		ln -sf "$DYNADJUST_INSTALL_PATH/libdnageoid.$_lib_ext"  "$BIN_FOLDER_FULLPATH/libdnageoid.$_lib_ext"
		sudo ln -sf "$DYNADJUST_INSTALL_PATH/libdnageoid.$_lib_ext" /opt/dynadjust/libdnageoid.$_lib_ext
		sudo ln -sf "$DYNADJUST_INSTALL_PATH/dnageoid" /opt/dynadjust/dnageoid
		echo " - dnageoid, libdnageoid.$_lib_ext"
	fi

	if [[ ( -e "../../bin/dnasegment" ) && ( "all" == "$_binary" || "$_binary" =~ "segment" ) ]]; then
		sudo cp ../../bin/libdnasegment.$_lib_ext "$DYNADJUST_INSTALL_PATH/"
		sudo cp ../../bin/dnasegment "$DYNADJUST_INSTALL_PATH/"
		ln -sf "$DYNADJUST_INSTALL_PATH/dnasegment" "$BIN_FOLDER_FULLPATH/dnasegment"
		ln -sf "$DYNADJUST_INSTALL_PATH/libdnasegment.$_lib_ext"  "$BIN_FOLDER_FULLPATH/libdnasegment.$_lib_ext"
		sudo ln -sf "$DYNADJUST_INSTALL_PATH/libdnasegment.$_lib_ext" /opt/dynadjust/libdnasegment.$_lib_ext 
		sudo ln -sf "$DYNADJUST_INSTALL_PATH/dnasegment" /opt/dynadjust/dnasegment
		echo " - dnasegment, libdnasegment.$_lib_ext"
	fi

	if [[ ( -e "../../bin/dnaplot" ) && ( "all" == "$_binary" || "$_binary" =~ "plot" ) ]]; then
		sudo cp ../../bin/libdnaplot.$_lib_ext "$DYNADJUST_INSTALL_PATH/"
		sudo cp ../../bin/dnaplot "$DYNADJUST_INSTALL_PATH/"
		ln -sf "$DYNADJUST_INSTALL_PATH/dnaplot" "$BIN_FOLDER_FULLPATH/dnaplot"
		ln -sf "$DYNADJUST_INSTALL_PATH/libdnaplot.$_lib_ext" "$BIN_FOLDER_FULLPATH/libdnaplot.$_lib_ext"
		sudo ln -sf "$DYNADJUST_INSTALL_PATH/libdnaplot.$_lib_ext" /opt/dynadjust/libdnaplot.$_lib_ext 
		sudo ln -sf "$DYNADJUST_INSTALL_PATH/dnaplot" /opt/dynadjust/dnaplot
		echo " - dnaplot, libdnaplot.$_lib_ext"
	fi
fi

# return to the original "current directory"
cd "$_cwd"

echo -e "Done.\n"

if [[ $_install -eq 0 ]]; then
	echo "Don't forget to add the bin directory to path in ~/.bash_profile"
	echo "For example:"
	echo -e "    EXPORT PATH=$PATH:$HOME/bin\n"
else
	echo "DynAdjust binaries can be located as follows:"
	if [[ -e "${_build_dir}/dynadjust/dynadjust/dynadjust" ]]; then
		echo " +[dynadjust]"
		echo "  ./dynadjust/build_gcc/dynadjust/dynadjust"
	fi
	if [[ -e "${_build_dir}/dynadjust/dnaadjustwrapper/dnaadjust" ]]; then
		echo " +[dnaadjust]"
		echo "  ./dynadjust/build_gcc/dynadjust/dnaadjustwrapper/dnaadjust"
		echo "  ./dynadjust/build_gcc/dynadjust/dnaadjust/libdnaadjust.$_lib_ext"
	fi

	if [[ -e "${_build_dir}/dynadjust/dnaimportwrapper/dnaimport" ]]; then
		echo " +[dnaimport]"
		echo "  ./dynadjust/build_gcc/dynadjust/dnaimport/libdnaimport.$_lib_ext"
		echo "  ./dynadjust/build_gcc/dynadjust/dnaimportwrapper/dnaimport"
	fi

	if [[ -e "${_build_dir}/dynadjust/dnareftranwrapper/dnareftran" ]]; then
		echo " +[dnareftran]"
		echo "  ./dynadjust/build_gcc/dynadjust/dnareftranwrapper/dnareftran"
		echo "  ./dynadjust/build_gcc/dynadjust/dnareftran/libdnareftran.$_lib_ext"
	fi

	if [[ -e "${_build_dir}/dynadjust/dnageoidwrapper/dnageoid" ]]; then
		echo " +[dnageoid]"
		echo "  ./dynadjust/build_gcc/dynadjust/dnageoidwrapper/dnageoid"
		echo "  ./dynadjust/build_gcc/dynadjust/dnageoid/libdnageoid.$_lib_ext"
	fi

	if [[ -e "${_build_dir}/dynadjust/dnasegmentwrapper/dnasegment" ]]; then
		echo " +[dnasegment]"
		echo "  ./dynadjust/build_gcc/dynadjust/dnasegmentwrapper/dnasegment"
		echo "  ./dynadjust/build_gcc/dynadjust/dnasegment/libdnasegment.$_lib_ext"
	fi

	if [[ -e "${_build_dir}/dynadjust/dnaplotwrapper/dnaplot" ]]; then
		echo " +[dnaplot]"
		echo "  ./dynadjust/build_gcc/dynadjust/dnaplotwrapper/dnaplot"
		echo "  ./dynadjust/build_gcc/dynadjust/dnaplot/libdnaplot.$_lib_ext"
	fi
	
	echo " "
fi
