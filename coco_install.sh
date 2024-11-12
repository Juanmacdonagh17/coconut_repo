#!/bin/bash

spinner() {
    # spinner function (just for fun)
    local msg="$1"
    shift
    local cmd=("$@")
    local pid
    "${cmd[@]}" &
    pid=$!

    local delay=0.1
    local spinstr='|/-\'
    local temp

    echo -n "$msg "

    while kill -0 "$pid" 2>/dev/null; do
        temp=${spinstr#?}
        printf "   %c  " "$spinstr"
        spinstr=$temp${spinstr%"$temp"}
        sleep $delay
        printf "\b\b\b\b\b\b"
    done

    wait "$pid"
    local exit_code=$?
    printf "    \b\b\b\b"

    if [ $exit_code -eq 0 ]; then
        echo " [DONE]"
    else
        echo " [FAILED]"
        exit $exit_code
    fi
}

echo 'checking for tools, curl.h and compiling coconut:'

# check for required tools
echo 'checking for gcc, make, and wget'

REQUIRED_TOOLS=(gcc wget make)
for tool in "${REQUIRED_TOOLS[@]}"; do
    if ! [ -x "$(command -v "$tool")" ]; then
        echo "Error: $tool is not installed. Install it and try again." >&2
        exit 1
    fi
done

# repo directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COCONUT_SOURCE="$SCRIPT_DIR/coconut.c"

# check for libcurl
if ! ldconfig -p | grep -q "libcurl\."; then
    echo "Library libcurl is not installed. Do you want to install it using the latest verison of Miniconda? (BE C A R E F U L) [y/n]"
    read -r response
    if [[ "$response" == [yY] ]]; then
        INSTALL_LIB=true
    else
        echo "Cannot proceed without libcurl. Exiting."
        exit 1
    fi
else
    # libcurl is already installed; get the path
    lib_paths=$(ldconfig -p | grep "libcurl\." | awk '{print $NF}' | sort | uniq)
    echo "libcurl is already installed at the following path(s):"
    echo "$lib_paths"
    INSTALL_LIB=false
    echo "compiling should go smoothly, if not try with miniconda :)"
fi

if [ "$INSTALL_LIB" = true ]; then
    echo 'installing Miniconda and setting up environment'

    # Install Miniconda
    mkdir -p "$SCRIPT_DIR/miniconda3"
    spinner "Downloading Miniconda" wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O "$SCRIPT_DIR/miniconda3/miniconda.sh"
    spinner "Installing Miniconda" bash "$SCRIPT_DIR/miniconda3/miniconda.sh" -b -u -p "$SCRIPT_DIR/miniconda3" > /dev/null 2>&1
    rm "$SCRIPT_DIR/miniconda3/miniconda.sh"

    # Initialize Conda
    export PATH="$SCRIPT_DIR/miniconda3/bin:$PATH"
    source "$SCRIPT_DIR/miniconda3/etc/profile.d/conda.sh"

    # create and activate conda env
    spinner "creating conda environment" conda create -y -n coconut_env
    conda activate coconut_env

    # install libcurl and gcc (gcc maybe it's already installed, but we do it anyways here)
    spinner "installing dependencies" conda install -y libcurl gcc_linux-64 gxx_linux-64

    # Set environment variables
   # CONDA_PREFIX=$(conda info --base)/envs/coconut_env
    export CPATH="$CONDA_PREFIX/include:$CPATH"
    export LIBRARY_PATH="$CONDA_PREFIX/lib:$LIBRARY_PATH"
    export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
fi

# compile coconut
echo 'compiling coconut:'

# set include and library paths if Miniconda was used
if [ "$INSTALL_LIB" = true ]; then
    # use gcc from the conda environment, not the default one
    GCC="$CONDA_PREFIX/bin/gcc"
    spinner "Compiling" "$GCC" -o coconut "$COCONUT_SOURCE" \
        -I"$CONDA_PREFIX/include" \
        -L"$CONDA_PREFIX/lib" \
        -Wl,-rpath,"$CONDA_PREFIX/lib" \
        -lcurl
else
    spinner "Compiling" gcc -o coconut "$COCONUT_SOURCE" -lcurl
fi

echo 'Installation complete :) To run coconut: ./coconut -help'

