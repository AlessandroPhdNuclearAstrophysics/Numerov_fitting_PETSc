#!/bin/bash

# Array to store missing packages
missing_packages=()

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

# Function to add missing packages to the list
add_missing_package() {
  if ! command -v "$1" &> /dev/null && ! check_package "$2"; then
    missing_packages+=("$2")
  fi
}

# Check for make
add_missing_package "make" "make"
if [[ "$OSTYPE" == "linux-gnu"* && ! dpkg -l | grep -q "make" ]]; then
  missing_packages+=("make")
fi

# Check for gfortran
add_missing_package "gfortran" "gfortran"

# Check for build-essential (Linux only)
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  add_missing_package "" "build-essential"
fi

# Check for lapack
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  add_missing_package "" "liblapack-dev"
elif [[ "$OSTYPE" == "darwin"* ]]; then
  if ! brew list lapack &> /dev/null; then
    missing_packages+=("lapack")
  fi
fi

# Check for blas
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  add_missing_package "" "blas"
elif [[ "$OSTYPE" == "darwin"* ]]; then
  if ! brew list openblas &> /dev/null; then
    missing_packages+=("openblas")
  fi
fi

# Check for xmgrace
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  add_missing_package "" "grace"
elif [[ "$OSTYPE" == "darwin"* ]]; then
  if ! brew list grace &> /dev/null; then
    missing_packages+=("grace")
  fi
fi

# Prompt for installation of all missing packages
if [[ ${#missing_packages[@]} -gt 0 ]]; then
  echo "The following packages are missing: ${missing_packages[*]}"
  read -p "Do you want to install them? (y/n): " choice
  if [[ "$choice" == "y" || "$choice" == "Y" ]]; then
    for package in "${missing_packages[@]}"; do
      case "$OSTYPE" in
        linux*)
          sudo apt-get update
          sudo apt-get install -y "$package"
          ;;
        darwin*)
          brew install "$package"
          ;;
        msys*|cygwin*|win32*)
          # Assuming Windows Subsystem for Linux (WSL) or compatible shell
          sudo apt-get update
          sudo apt-get install -y "$package"
          ;;
        *)
          echo "Unsupported OS: $OSTYPE"
          exit 1
          ;;
      esac
    done
  else
    echo "Skipping installation of missing packages."
  fi
else
  echo "All required packages are installed."
fi
