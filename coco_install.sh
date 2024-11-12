#!/bin/bash

spinner() {
    # spinner 
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

echo 'Installing coconut and curl.h'

# check for required tools
echo 'checking for gcc, make and wget'

REQUIRED_TOOLS=(gcc wget make)
for tool in "${REQUIRED_TOOLS[@]}"; do
    if ! [ -x "$(command -v "$tool")" ]; then
        echo "$tool is not installed. Install it and try again :)" >&2
        exit 1
    fi
done

# dir
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COCONUT_SOURCE="$SCRIPT_DIR/coconut.c"

# install dependencies locally
mkdir -p "$HOME/local" # maybe make this inside the repo ? so no to make an extra local folder? idk
cd "$HOME/local" || exit 

# download and install libcurl
echo 'downloading and installing libcurl \n'

# downloading curl
spinner "downloading curl" wget https://curl.se/download/curl-7.88.1.tar.gz

# extracting curl
spinner "extracting curl" tar -xzf curl-7.88.1.tar.gz

cd curl-7.88.1 || exit

# configuring curl
spinner "configuring curl" ./configure --prefix="$HOME/local" > /dev/null 2>&1

# compiling curl
spinner "compiling curl" make > /dev/null 2>&1

# installing curl
spinner "installing curl" make install > /dev/null 2>&1

# compile coconut
spinner "compiling coconut" gcc -o coconut "$COCONUT_SOURCE" -I"$HOME/local/include" -L"$HOME/local/lib" -lcurl

echo 'installation complete :D'

