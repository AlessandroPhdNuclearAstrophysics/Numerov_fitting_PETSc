# Project: Fit Numerov with PETSc

The program performs the following tasks:

1. **Quadratic Interpolation**: Computes a quadratic interpolation of the observable $k^3 \cot\delta$ as a function of $k^2$ for center-of-mass energies ($E_{\text{cm}}$) ranging from 0.05 to 1 MeV, using the AV18 potential in the 1P1 channel.

2. **Parameter Fitting**: Fits the depth ($C$) and cutoff ($R_0$) of the EFT-pless potential to match the quadratic interpolation parameters ($a$, $b$, and $c$) derived from the AV18 potential, where $k^3 \cot \delta = a k^4 + b k^2 + c$.

3. **Scattering Properties**: Evaluates the scattering volume and the effective range based on the fitted EFT-pless potential.

This process enables a detailed comparison between the AV18 and EFT-pless potentials and provides insights into the scattering behavior in the 1P1 channel.

## Project Structure
The project directory is organized as follows:
- **`bin/`**: Directory containing the compiled executable file.
- **`build/`**: Directory containing object files (`.o`) generated during the compilation process.
- **`src/`**: Source code files implementing the Numerov algorithm and data fitting.
- **`output/`**: Directory for generated output files, including numerical results and plots.
- **`Makefile`**: Build system to compile, execute, and visualize results.

## Output Files
The `output/` directory contains the following files:

- **`compare.agr`**: A file in the `AGR` format for graphical comparison of results.
- **`kcotd_1P1.dat`**: Numerical data file containing computed results for the observable $k^3 \cot\delta$ for the channel 1P1, generated using the EFT-pless fitted potential.
- **`kcotd_1P1_data.dat`**: Numerical data file containing computed results for the observable $k^3 \cot\delta$ for the channel 1P1, generated using the AV18 potential.

These files are generated during the execution of the program and can be used for further analysis or visualization.

## Usage Instructions

### 1. Compile the Project
To compile the source code, run:
```bash
make
```

### 2. Execute the Program
Run the compiled program with:
```bash
make run
```
This command runs the executable finding the fitting parameters and generates the output files in the `output/` directory.

### 3. Visualize the Results
To display the graphical plot of the results, use:
```bash
make compare
```
This will open the generated plot (`output/compare.agr`).

### 4. Clean the Project
To remove compiled files, execute:
```bash
make clean
```

### 5. Delete Output Data Files
To remove all `.dat` files from the `output/` directory, use:
```bash
make del_out
```
This command deletes the numerical data files generated during program execution.

## Dependencies
Ensure the following dependencies are installed on your system:
- PETSc library
- GNU Make
- XmGrace (required for generating plots)

## Additional Notes
- The `Makefile` can be modified to adjust paths, compiler options, or dependencies as needed.
- For further details, refer to the comments in the source code files located in the `src/` directory.
- Ensure all dependencies are correctly configured to avoid runtime errors.

