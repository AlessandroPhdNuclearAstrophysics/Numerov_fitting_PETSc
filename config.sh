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
missing_packages_list=()

# Function to check if a package is installed
check_package() {
  case "$OSTYPE" in
    linux*)
      dpkg -l | grep -q "$1"
      ;;
    darwin*)
      brew list "$1" &> /dev/null
      ;;
    *)
      echo "Unsupported OS: $OSTYPE"
      exit 1
      ;;
  esac
}

# Function to add missing package flag
add_missing_package_flag() {
  local command_name="$1"
  local package_name="$2"
  local flag="$3"

  if ! command -v "$command_name" &> /dev/null && ! check_package "$package_name"; then
    missing_packages_binary=$((missing_packages_binary | flag))
    missing_packages_list+=("$package_name")
  fi
}

# Function to install a package
install_package() {
  local package="$1"
  local flag="$2"

  case "$OSTYPE" in
    linux*)
      if sudo apt-get install -y "$package"; then
        missing_packages_binary=$((missing_packages_binary & ~flag))
      fi
      ;;
    darwin*)
      if brew install "$package"; then
        missing_packages_binary=$((missing_packages_binary & ~flag))
      fi
      ;;
    *)
      echo "Unsupported OS: $OSTYPE"
      exit 1
      ;;
  esac
}

# Check for required packages
add_missing_package_flag "make" "make" $MAKE_FLAG
add_missing_package_flag "gfortran" "gfortran" $GFORTRAN_FLAG

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  add_missing_package_flag "" "build-essential" $BUILD_ESSENTIAL_FLAG
  add_missing_package_flag "" "liblapack-dev" $LAPACK_FLAG
  add_missing_package_flag "" "blas" $BLAS_FLAG
  add_missing_package_flag "" "grace" $GRACE_FLAG
elif [[ "$OSTYPE" == "darwin"* ]]; then
  add_missing_package_flag "" "lapack" $LAPACK_FLAG
  add_missing_package_flag "" "openblas" $BLAS_FLAG
  add_missing_package_flag "" "grace" $GRACE_FLAG
fi

echo "Missing packages binary: $missing_packages_binary"

# Prompt user to install missing packages
if [[ ${#missing_packages_list[@]} -gt 0 ]]; then
  echo "The following packages are missing: ${missing_packages_list[*]}"
  read -p "Do you want to install all of them? (y/n): " response
  if [[ "$response" =~ ^[yY]$ ]]; then
    for package in "${missing_packages_list[@]}"; do
      case "$package" in
        "make") install_package "make" $MAKE_FLAG ;;
        "gfortran") install_package "gfortran" $GFORTRAN_FLAG ;;
        "build-essential") install_package "build-essential" $BUILD_ESSENTIAL_FLAG ;;
        "liblapack-dev") install_package "liblapack-dev" $LAPACK_FLAG ;;
        "blas") install_package "blas" $BLAS_FLAG ;;
        "grace") install_package "grace" $GRACE_FLAG ;;
        "lapack") install_package "lapack" $LAPACK_FLAG ;;
        "openblas") install_package "openblas" $BLAS_FLAG ;;
      esac
    done
  else
    echo "You can choose specific packages to install."
    for i in "${!missing_packages_list[@]}"; do
      echo "$((i + 1))) ${missing_packages_list[i]}"
    done
    read -p "Enter the numbers of the packages you want to install (space-separated): " selected_numbers
    for number in $selected_numbers; do
      if [[ $number -ge 1 && $number -le ${#missing_packages_list[@]} ]]; then
        package="${missing_packages_list[$((number - 1))]}"
        case "$package" in
          "make") install_package "make" $MAKE_FLAG ;;
          "gfortran") install_package "gfortran" $GFORTRAN_FLAG ;;
          "build-essential") install_package "build-essential" $BUILD_ESSENTIAL_FLAG ;;
          "liblapack-dev") install_package "liblapack-dev" $LAPACK_FLAG ;;
          "blas") install_package "blas" $BLAS_FLAG ;;
          "grace") install_package "grace" $GRACE_FLAG ;;
          "lapack") install_package "lapack" $LAPACK_FLAG ;;
          "openblas") install_package "openblas" $BLAS_FLAG ;;
        esac
      else
        echo "Invalid selection: $number"
      fi
    done
  fi
fi

echo "Missing packages binary after installation: $missing_packages_binary"

# Final check
if [[ $missing_packages_binary -eq 0 ]]; then
  echo "All required packages are installed."
else
  echo "Some required packages are still missing."
fi

exit $missing_packages_binary
