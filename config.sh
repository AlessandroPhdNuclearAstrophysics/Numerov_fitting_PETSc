#!/bin/bash

# Binary flags for each package
MAKE_FLAG=1          # 000001
GFORTRAN_FLAG=2      # 000010
BUILD_ESSENTIAL_FLAG=4 # 000100
LAPACK_FLAG=8        # 001000
BLAS_FLAG=16         # 010000
GRACE_FLAG=32        # 100000

# Variable to store missing packages as a binary value
missing_packages_binary=0

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
      dpkg -l | grep -q "$1"
      ;;
    *)
      echo "Unsupported OS: $OSTYPE"
      exit 1
      ;;
  esac
}

# Function to add missing package flag
add_missing_package_flag() {
  if ! command -v "$1" &> /dev/null && ! check_package "$2"; then
    missing_packages_binary=$((missing_packages_binary | $3))
  fi
}

# Check for make
add_missing_package_flag "make" "make" $MAKE_FLAG

# Check for gfortran
add_missing_package_flag "gfortran" "gfortran" $GFORTRAN_FLAG

# Check for build-essential (Linux only)
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  add_missing_package_flag "" "build-essential" $BUILD_ESSENTIAL_FLAG
fi

# Check for lapack
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  add_missing_package_flag "" "liblapack-dev" $LAPACK_FLAG
elif [[ "$OSTYPE" == "darwin"* ]]; then
  if ! brew list lapack &> /dev/null; then
    missing_packages_binary=$((missing_packages_binary | $LAPACK_FLAG))
  fi
fi

# Check for blas
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  add_missing_package_flag "" "blas" $BLAS_FLAG
elif [[ "$OSTYPE" == "darwin"* ]]; then
  if ! brew list openblas &> /dev/null; then
    missing_packages_binary=$((missing_packages_binary | $BLAS_FLAG))
  fi
fi

# Check for xmgrace
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  add_missing_package_flag "" "grace" $GRACE_FLAG
elif [[ "$OSTYPE" == "darwin"* ]]; then
  if ! brew list grace &> /dev/null; then
    missing_packages_binary=$((missing_packages_binary | $GRACE_FLAG))
  fi
fi

# Return the binary flag for missing packages
exit $missing_packages_binary
