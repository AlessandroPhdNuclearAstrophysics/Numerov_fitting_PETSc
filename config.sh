#!/bin/bash

# Function to check if a package is installed
check_package() {
  case "$OSTYPE" in
    linux*)
      dpkg -l | grep -q "$1"
      ;;
    darwin*)
      brew list "$1" &> /dev/null
      ;;
    msys*|cygwin*|win32*)
      # Assuming Windows Subsystem for Linux (WSL) or compatible shell
      dpkg -l | grep -q "$1"
      ;;
    *)
      echo "Unsupported OS: $OSTYPE"
      exit 1
      ;;
  esac
}

# Function to prompt for installation
prompt_install() {
  read -p "$1 is not installed. Do you want to install it? (y/n): " choice
  if [[ "$choice" == "y" || "$choice" == "Y" ]]; then
    case "$OSTYPE" in
      linux*)
        sudo apt-get update
        sudo apt-get install -y "$1"
        ;;
      darwin*)
        brew install "$1"
        ;;
      msys*|cygwin*|win32*)
        # Assuming Windows Subsystem for Linux (WSL) or compatible shell
        sudo apt-get update
        sudo apt-get install -y "$1"
        ;;
      *)
        echo "Unsupported OS: $OSTYPE"
        exit 1
        ;;
    esac
  else
    echo "Skipping installation of $1."
  fi
}

# Check for make
if ! command -v make &> /dev/null; then
  prompt_install "make"
else
  echo "make is installed."
fi

# Check for gfortran
if ! command -v gfortran &> /dev/null; then
  prompt_install "gfortran"
else
  echo "gfortran is installed."
fi

# Check for build-essential (Linux only)
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  if ! check_package "build-essential"; then
    prompt_install "build-essential"
  else
    echo "build-essential is installed."
  fi
fi

# Check for lapack
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  if ! check_package "liblapack-dev"; then
    prompt_install "liblapack-dev"
  else
    echo "liblapack-dev (lapack) is installed."
  fi
elif [[ "$OSTYPE" == "darwin"* ]]; then
  if ! brew list lapack &> /dev/null; then
    prompt_install "lapack"
  else
    echo "lapack is installed."
  fi
fi

# Check for blas
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  if ! check_package "blas"; then
    prompt_install "blas"
  else
    echo "blas is installed."
  fi
elif [[ "$OSTYPE" == "darwin"* ]]; then
  if ! brew list openblas &> /dev/null; then
    prompt_install "openblas"
  else
    echo "openblas (blas) is installed."
  fi
fi

# Check for xmgrace
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  if ! check_package "grace"; then
    prompt_install "grace"
  else
    echo "xmgrace (grace) is installed."
  fi
elif [[ "$OSTYPE" == "darwin"* ]]; then
  if ! brew list grace &> /dev/null; then
    prompt_install "grace"
  else
    echo "xmgrace (grace) is installed."
  fi
fi