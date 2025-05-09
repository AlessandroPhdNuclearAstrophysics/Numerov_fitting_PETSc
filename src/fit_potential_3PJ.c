/*
   Include "petsctao.h" so that we can use TAO solvers.  Note that this
   file automatically includes libraries such as:
     petsc.h       - base PETSc routines   petscvec.h - vectors
     petscsys.h    - system routines       petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
*/

#include <petsctao.h>
#include "libs/scattering_single_channel.h"
#include "libs/physical_constants.h"
#include "libs/read_input_file.h"

static char help[] = "Find the constants to fit the 3PJ channel evaluated \n\
            with AV18, but using the EFT-pless model \n";

#define NOBSERVATIONS 6
#define NPARAMETERS   3
#define NCHANNELS NOBSERVATIONS/3


const int ipot_ref = 18;
const int ilb_ref  = 1;
const int ipot = -1;
const int ilb  = -1;
const int ne = 20;
const double emin = 0.0;
const double emax = 1.0;
const double he   = (emax-emin)/ne;
double *energies;
double *k2;
double *kcotd[NCHANNELS];

PetscInt J_ref[] = {0, 1};
PetscInt J_fit[] = {0, 1};

QuantumNumbers qn_ref = {18,  1, 1, 1, J_ref, 2};
QuantumNumbers qn_fit = {-1, -1, 1, 1, J_fit, 2};

LECs lecs = {0.0, NULL, 4};

DataRow data;



/* User-defined application context */
typedef struct {
  /* Working space */
  PetscReal t[NOBSERVATIONS];              /* array of independent variables of observation */
  PetscReal y[NOBSERVATIONS];              /* array of dependent variables */
  PetscReal j[NOBSERVATIONS][NPARAMETERS]; /* dense jacobian matrix array*/
  PetscInt  idm[NOBSERVATIONS];            /* Matrix indices for jacobian */
  PetscInt  idn[NPARAMETERS];
} AppCtx;

/* User provided Routines */
PetscErrorCode InitializeData(AppCtx *user);
PetscErrorCode FormStartingPoint(Vec);
PetscErrorCode EvaluateFunction(Tao, Vec, Vec, void *);
PetscErrorCode EvaluateJacobian(Tao, Vec, Mat, Mat, void *);
PetscErrorCode ComputeParameterUncertainty(Mat, Vec);


/*--------------------------------------------------------------------*/

/* argc and argv = standard parameters in a C program's main function.

- argc is the argument count, representing the number of command-line arguments 
(including the program name itself).
- argv is an array of strings (character pointers) containing these arguments.

They are passed to PetscInitialize to allow the PETSc library to parse and handle 
any PETSc-specific command-line arguments, such as options for solver configuration, 
logging, and debugging.
PETSc often provides command-line options like -tao_monitor or -tao_type for 
tuning the behavior of solvers without modifying the code.*/

/*--------------------------------------------------------------------*/
int main(int argc, char **argv)
{
  Vec       x, xu, xl, f; /* solution, function */
  Tao       tao;  /* Tao solver context */
  PetscReal hist[100], resid[100];
  PetscInt  lits[100];
  AppCtx    user; /* user-defined work context */
  

  /* Preparing the values of energies to evaluate k^3 cot(delta) and 
      therefore the quadratic fitting parameters a, b and c. */
  energies = (double *) malloc(sizeof(double)*ne);
  for (int i=0; i<NCHANNELS; i++) {
    kcotd[i] = (double *) malloc(sizeof(double)*ne);
  }
  for (int ie=1; ie<= ne; ie++) energies[ie-1] = he*ie;

  lecs.LECS = (double *) malloc(sizeof(double)*lecs.NLECS);
  lecs.LECS[0] = 0.0;
  lecs.LECS[1] = 0.0;
  lecs.LECS[2] = -0.136793827E1;
  lecs.LECS[3] =  0.219960910E-1;


  PetscFunctionBeginUser;
  /* MACRO PetscFunctionBeginUser:
    Functionality:
    - It initializes the function's scope for PETSc-related debugging or profiling.
    - Ensures that any errors encountered in the function are correctly logged and handled.
    Distinction from PetscFunctionBegin:
    - While PetscFunctionBegin is the general-purpose version used internally in PETSc,
     PetscFunctionBeginUser is meant specifically for user-defined routines, helping PETSc 
     differentiate between user functions and library functions. */
  PetscCall(PetscInitialize(&argc, &argv, NULL, help));

  if (read_file("input_values_3PJ.dat") != 0) {
    PetscPrintf(PETSC_COMM_SELF, "Error reading file\n");
    return 1;
  }
  PetscPrintf(PETSC_COMM_SELF, "Input file read file successfully\n");
  const int num_starting_points = get_total_rows();
  if (num_starting_points < 1) {
    PetscPrintf(PETSC_COMM_SELF, "No data found in the file\n");
    return 1;
  }
  PetscPrintf(PETSC_COMM_SELF, "Number of starting points: %d\n", num_starting_points);

   /* Allocate vectors */
  PetscCall(VecCreateSeq(MPI_COMM_SELF, NPARAMETERS, &x));
  PetscCall(VecCreateSeq(MPI_COMM_SELF, NOBSERVATIONS, &f));

  /* Create TAO solver and set desired solution method */
  PetscCall(TaoCreate(PETSC_COMM_SELF, &tao));
  PetscCall(TaoSetType(tao, TAOPOUNDERS));

  /* Creates the 2 dimensional vectors to set the upper and lower bounds 
      for the parameters */
  VecCreateSeq(PETSC_COMM_SELF, NPARAMETERS, &x);
  VecDuplicate(x, &xl);  // limiti inferiori
  VecDuplicate(x, &xu);  // limiti superiori

  /* Set the limits on R (cutoff) and C(depth of the potential) */
  PetscScalar lower[3] = {0.5, -20, -20};  // r min, c min
  PetscScalar upper[3] = {3.5,  20,  20};  // r max, c max

  VecSetValues(xl, 3, (PetscInt[]){0,1,2}, lower, INSERT_VALUES);
  VecSetValues(xu, 3, (PetscInt[]){0,1,2}, upper, INSERT_VALUES);
  VecAssemblyBegin(xl); VecAssemblyEnd(xl);
  VecAssemblyBegin(xu); VecAssemblyEnd(xu);

  /* Associates these limits to TAO */
  TaoSetVariableBounds(tao, xl, xu);

  for (int i_point=0; i_point<num_starting_points; i_point++) {
    /* Get the starting point from the file */
    data = get_row(i_point);

    /* Set the function routines. */
    PetscCall(InitializeData(&user));
    PetscCall(FormStartingPoint(x));
    PetscCall(TaoSetSolution(tao, x));
    PetscCall(TaoSetResidualRoutine(tao, f, EvaluateFunction, (void *)&user));

    /* Check for any TAO command line arguments */
    PetscCall(TaoSetFromOptions(tao));
    PetscCall(TaoSetConvergenceHistory(tao, hist, resid, 0, lits, 100, PETSC_TRUE));

    /* Perform the Solve */
    PetscCall(TaoSolve(tao));
    
    /* Recast the solution to an array to print and use the solution */
    const PetscScalar *x_array;
    VecGetArrayRead(x, &x_array);
    const double sol[] = { x_array[0], x_array[1], x_array[2] };
    VecRestoreArrayRead(x, &x_array);

    lecs.R = sol[0];
    lecs.LECS[0] = sol[1];
    lecs.LECS[1] = sol[2];
    Observables obs = scattering_numerov_quantum_num(energies, kcotd, ne, &qn_fit, &lecs);
    PetscPrintf(PETSC_COMM_SELF, "\n\n----------------------------------------\n");
    PetscPrintf(PETSC_COMM_SELF, "SOLUTION\n");
    PetscPrintf(PETSC_COMM_SELF, "----------------------------------------\n");
    PetscPrintf(PETSC_COMM_SELF, "R: %.15f C11: %.15f C6: %.15f\n", sol[0], sol[1], sol[2]);
    for (int i=0; i<NCHANNELS; i++){
      PetscPrintf(PETSC_COMM_SELF, "c: %.15f\tb: %.15f\ta: %.15f\n", obs.coeffs[i].cba[0],obs.coeffs[i].cba[1],obs.coeffs[i].cba[2]);
      PetscPrintf(PETSC_COMM_SELF, "Scattering length: %f15\n", -1/obs.coeffs[i].cba[0]);
      PetscPrintf(PETSC_COMM_SELF, "Effective range  : %f15\n\n\n", 2*obs.coeffs[i].cba[1]);
    }

    /* Print everything to file */
    FILE *fp  = fopen(data.output_filename,"write");
    if (fp==NULL) {
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Error opening file");
      return 1;
    }
    
    fprintf(fp,"#R: %.15f C11: %.15f C6: %.15f\n", sol[0], sol[1], sol[2]);
    for (int i=0; i<NCHANNELS; i++) {
      fprintf(fp,"#c: %.15f\tb: %.15f\ta: %.15f\n", obs.coeffs[i].cba[0],obs.coeffs[i].cba[1],obs.coeffs[i].cba[2]);
      fprintf(fp,"#Scattering length: %.15f\n", -1/obs.coeffs[0].cba[0]);
      fprintf(fp,"#Effective range  : %.15f\n\n\n", 2*obs.coeffs[0].cba[1]);
      for (int ie=1; ie <= ne; ie++) fprintf(fp, "%.15f\t%.15f\n", k2_from_E(energies[ie-1]), kcotd[i][ie-1]);
    }

    fclose(fp);
  }

  /* Free TAO data structures */
  PetscCall(TaoDestroy(&tao));

  /* Free PETSc data structures */
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&f));
  PetscCall(VecDestroy(&xl));
  PetscCall(VecDestroy(&xu));

  PetscCall(PetscFinalize());
  return 0;
} 























/*--------------------------------------------------------------------*/
PetscErrorCode EvaluateFunction(Tao tao, Vec X, Vec F, void *ptr)
{
  AppCtx          *user = (AppCtx *)ptr; // Cast the context pointer to the user-defined context
  PetscInt         i; // Loop index
  const PetscReal *x; // Pointer to the array of parameters
  PetscReal       *y = user->y, *f; // Pointers to the arrays of observations and independent variables

  PetscFunctionBegin; // Begin PETSc function
  PetscCall(VecGetArrayRead(X, &x)); // Get read-only access to the array inside vector X
  PetscCall(VecGetArray(F, &f)); // Get access to the array inside vector F

  double sum_f = 0.0;
  // Compute the residuals for each observation
  lecs.R = x[0];
  lecs.LECS[0] = x[1];
  lecs.LECS[1] = x[2];
  Observables obs = scattering_numerov_quantum_num(energies, kcotd, ne, &qn_fit, &lecs);
  double weight[] = { 5, 3, 0};
  
  PetscPrintf(PETSC_COMM_SELF, "\n\nR: %.15f C11: %.15f C6: %.15f\n", x[0], x[1], x[2]);
  for (int ich=0; ich < NCHANNELS; ich++) {
    coeffs coeff = obs.coeffs[ich];
    for (i = 0; i < 3; i++) {
      f[ich*3 + i] = (y[ich*3 + i] - coeff.cba[i])/(y[ich*3 + i] + coeff.cba[i]);
      f[ich*3 + i] *= weight[i];
    }

    for (i=0; i < NOBSERVATIONS; i++) sum_f += f[i]*f[i];
    PetscPrintf(PETSC_COMM_SELF, "c: %.15f\tb: %.15f\ta: %.15f\n", coeff.cba[0],coeff.cba[1],coeff.cba[2]);
  }
  PetscPrintf(PETSC_COMM_SELF, "sum_f: %.15f\n", sum_f);

  


  PetscCall(VecRestoreArrayRead(X, &x)); // Restore the array inside vector X
  PetscCall(VecRestoreArray(F, &f)); // Restore the array inside vector F
// The function PetscLogFlops is used to log the number of floating-point operations (FLOPs) 
// performed by the program. This is useful for performance monitoring and optimization, as it 
//allows developers to track the computational workload of their code.
  PetscCall(PetscLogFlops(6 * NOBSERVATIONS)); // Log the number of floating-point operations
  PetscFunctionReturn(PETSC_SUCCESS); // Return success
}

/*------------------------------------------------------------*/
/* J[i][j] = df[i]/dt[j] */
PetscErrorCode EvaluateJacobian(Tao tao, Vec X, Mat J, Mat Jpre, void *ptr)
{
  AppCtx          *user = (AppCtx *)ptr;
  PetscInt         i;
  const PetscReal *x;
  PetscReal       *t = user->t;
  PetscReal        base;

  PetscFunctionBegin;
//  PetscPrintf(PETSC_COMM_SELF, "EvaluateJacobian\n");
  PetscCall(VecGetArrayRead(X, &x));
  for (i = 0; i < NOBSERVATIONS; i++) {
    base = PetscExpScalar(-x[0] * t[i]) / (x[1] + x[2] * t[i]);
//    PetscPrintf(PETSC_COMM_SELF, "base: %f\n", base);

    user->j[i][0] = t[i] * base;
    user->j[i][1] = base / (x[1] + x[2] * t[i]);
    user->j[i][2] = base * t[i] / (x[1] + x[2] * t[i]);
  }

  /* Assemble the matrix */
  PetscCall(MatSetValues(J, NOBSERVATIONS, user->idm, NPARAMETERS, user->idn, (PetscReal *)user->j, INSERT_VALUES));
  PetscCall(MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY));

  PetscCall(VecRestoreArrayRead(X, &x));
  PetscCall(PetscLogFlops(NOBSERVATIONS * 13));
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* ------------------------------------------------------------ */
PetscErrorCode FormStartingPoint(Vec X)
{
  PetscReal *x;

  PetscFunctionBegin;
  PetscCall(VecGetArray(X, &x));
  x[0] =  data.R;    
  x[1] =  data.C00;
  x[2] =  data.C6;
  lecs.R = x[0];
  lecs.LECS[0] = x[1];
  lecs.LECS[1] = x[2];
//   x[0] =  0.5; // Try these to see that changing to a different number failes to make it converge
//   x[1] = -1;   // Try these to see that changing to a different number failes to make it converge
// x[0] = 0.7;    // Try these to see that changing even a little failes to find the best point
// x[1] = -2.5;   // Try these to see that changing even a little failes to find the best point
  PetscCall(VecRestoreArray(X, &x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* ------------------------------------------------------------ */
/* Function to compute and print parameter uncertainties - aggiunta da me*/
PetscErrorCode ComputeParameterUncertainty(Mat J, Vec f)
{
  Mat            JTJ;     // J^T * J matrix
  Vec            diag;    // Diagonal vector for variances
  PetscReal      norm_f;  // Norm of residual vector
  PetscReal      sigma2;  // Variance of residuals
  PetscReal     *uncert;  // Array to store uncertainties
  KSP            ksp;     // Krylov Subspace Solver
  PC             pc;      // Preconditioner context
  Vec            e_i, cov_i; // Unit vector and covariance column

  PetscFunctionBegin;
  // Compute J^T * J
  PetscCall(MatTransposeMatMult(J, J, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &JTJ));

  // Compute norm of residual vector f
  PetscCall(VecNorm(f, NORM_2, &norm_f));
  sigma2 = (norm_f * norm_f) / (NOBSERVATIONS - NPARAMETERS);

  // Create a vector for the diagonal of the covariance matrix
  PetscCall(VecCreateSeq(PETSC_COMM_SELF, NPARAMETERS, &diag));
  PetscCall(VecSet(diag, 0.0)); // Initialize to zero

  // Set up KSP solver to solve JTJ * x = b
  PetscCall(KSPCreate(PETSC_COMM_SELF, &ksp));
  PetscCall(KSPSetOperators(ksp, JTJ, JTJ));
  PetscCall(KSPGetPC(ksp, &pc));
  PetscCall(PCSetType(pc, PCLU)); // Use LU factorization
  PetscCall(KSPSetFromOptions(ksp));

  // Solve for each diagonal entry
  PetscCall(VecCreateSeq(PETSC_COMM_SELF, NPARAMETERS, &e_i));
  PetscCall(VecCreateSeq(PETSC_COMM_SELF, NPARAMETERS, &cov_i));
  for (PetscInt i = 0; i < NPARAMETERS; i++) {
    // Set e_i to the i-th unit vector
    PetscCall(VecSet(e_i, 0.0));
    PetscCall(VecSetValue(e_i, i, 1.0, INSERT_VALUES));
    PetscCall(VecAssemblyBegin(e_i));
    PetscCall(VecAssemblyEnd(e_i));

    // Solve JTJ * cov_i = e_i
    PetscCall(KSPSolve(ksp, e_i, cov_i));

    // Extract the i-th diagonal element (variance)
    PetscReal var_i;
    PetscCall(VecGetValues(cov_i, 1, &i, &var_i));
    PetscCall(VecSetValue(diag, i, var_i, INSERT_VALUES));
  }
  PetscCall(VecAssemblyBegin(diag));
  PetscCall(VecAssemblyEnd(diag));

  // Compute uncertainties (square roots of variances)
  PetscCall(VecGetArray(diag, &uncert));
  for (PetscInt i = 0; i < NPARAMETERS; i++) {
    uncert[i] = PetscSqrtReal(sigma2 * uncert[i]);
    PetscPrintf(PETSC_COMM_SELF, "Uncertainty in parameter %d: %g\n", i, uncert[i]);
  }
  PetscCall(VecRestoreArray(diag, &uncert));

  // Clean up
  PetscCall(VecDestroy(&e_i));
  PetscCall(VecDestroy(&cov_i));
  PetscCall(VecDestroy(&diag));
  PetscCall(MatDestroy(&JTJ));
  PetscCall(KSPDestroy(&ksp));

  PetscFunctionReturn(PETSC_SUCCESS);
}



/* ---------------------------------------------------------------------- */
PetscErrorCode InitializeData(AppCtx *user)
{
  PetscReal *y = user->y;

  PetscPrintf(PETSC_COMM_SELF, "\n\n\nEvaluating data to be fitted with AV18...\n");
  Observables obs = scattering_numerov_quantum_num(energies, kcotd, ne, &qn_ref, &lecs);

  for (int ich=0; ich < NCHANNELS; ich++) {
    coeffs coeff = obs.coeffs[ich];
    PetscPrintf(PETSC_COMM_SELF, "x[0]: %.15f x[1]: %.15f x[2]: %.15f\n", coeff.cba[0], coeff.cba[1], coeff.cba[2]);
    PetscPrintf(PETSC_COMM_SELF, "#Scattering length: %.15f\n", -1/coeff.cba[0]);
    PetscPrintf(PETSC_COMM_SELF, "#Effective range  : %.15f\n\n\n", 2*coeff.cba[1]);
    y[ich * 3 + 0] = coeff.cba[0];
    y[ich * 3 + 1] = coeff.cba[1];
    y[ich * 3 + 2] = coeff.cba[2];
  }
  FILE *fp  = fopen("output/kcotd_3PJ_data.dat","write");
  if (fp==NULL) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Error opening file");
    return 1;
  }
  fprintf(fp,"# --- potential data to be fitted are a, b and c ---\n");
  for (int ich=0; ich < NCHANNELS; ich++) {
    coeffs coeff = obs.coeffs[ich];
    fprintf(fp,"# --- channel %d ---\n", ich);
    fprintf(fp,"#c: %.15f\tb: %.15f\ta: %.15f\n", coeff.cba[0],coeff.cba[1],coeff.cba[2]);
    fprintf(fp,"#Scattering length: %.15f\n", -1/coeff.cba[0]);
    fprintf(fp,"#Effective range  : %.15f\n\n\n", 2*coeff.cba[1]);
    for (int ie=1; ie <= ne; ie++) fprintf(fp, "%.15f\t%.15f\n", k2_from_E(energies[ie-1]), kcotd[ich][ie-1]);
  }
  fclose(fp);
  PetscFunctionBegin;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   build:
      requires: !complex !single

1)   test:
      args: -tao_monitor_short -tao_max_it 100 -tao_type pounders -tao_gatol 1.e-5
      -tao_view

features:
-tao_monitor_short: Enables short monitoring of the TAO solver's progress.
-tao_max_it 100: Sets the maximum number of iterations to 100.
-tao_type pounders: Uses the POUNDERS (Practical Optimization Using No Derivatives for 
sums of Squares) solver type. It is a derivative-free optimization method implemented in 
PETSc's TAO (Toolkit for Advanced Optimization).
Use Case: POUNDERS is suitable for problems where the objective function is a sum of squares 
and derivatives are unavailable or expensive to compute.
Key Feature: It employs quadratic models to approximate the objective function and iteratively 
minimizes it, typically used for black-box optimization.
-tao_gatol 1.e-5: Sets the gradient tolerance for convergence to 1.e-5.
-tao_view: Displays the TAO solver's configuration and results.

2)   test:
      suffix: 2
      args: -tao_monitor_short -tao_max_it 100 -tao_type brgn 
      -tao_brgn_regularization_type l2prox -tao_brgn_regularizer_weight 1e-4 -tao_gatol 1.e-5

features:
-tao_monitor_short: Enables short monitoring of the TAO solver's progress.
-tao_max_it 100: Sets the maximum number of iterations to 100.
-tao_type brgn: Uses the BRGN (Bounded Regularized Gauss-Newton) solver type.
BRGN stands for Bounded Regularized Gauss-Newton.
This method is designed for nonlinear least-squares problems and incorporates regularization 
to stabilize the solution. L2 Prox refers to L2 -regularization (proximal term), which adds 
a quadratic penalty on the parameters to prevent large variations and ensure smoother solutions.
Key Features:
Uses Gauss-Newton approximations for solving nonlinear least squares problems.
The L2-proximal term helps in controlling the magnitude of the variables.
-tao_brgn_regularization_type l2prox: Uses L2 proximal regularization.
-tao_brgn_regularizer_weight 1e-4: Sets the regularization weight to 1e-4.
-tao_gatol 1.e-5: Sets the gradient tolerance for convergence to 1.e-5.

3)   test:
      suffix: 3
      args: -tao_monitor_short -tao_max_it 100 -tao_type brgn -tao_brgn_regularization_type 
      l1dict -tao_brgn_regularizer_weight 1e-4 -tao_brgn_l1_smooth_epsilon 1e-6 -tao_gatol 
      1.e-5

features:
-tao_monitor_short: Enables short monitoring of the TAO solver's progress.
-tao_max_it 100: Sets the maximum number of iterations to 100.
-tao_type brgn: Uses the BRGN (Bounded Regularized Gauss-Newton) solver type.
This variant of BRGN uses L1-regularization, promoting sparsity in the solution.
The term Dict likely refers to the use of a dictionary approach, where the solution is 
expressed as a sparse combination of basis elements.
Use Case: Useful for problems requiring sparse solutions (e.g., feature selection or 
compressed sensing). Vuol dire che la soluzione è un vettore con molte entrate = 0.
Key Features: Encourages sparsity in the optimization process via L1-norm.
Often used when the underlying solution has many zero or near-zero components.
-tao_brgn_regularization_type l1dict: Uses L1 dictionary regularization.
-tao_brgn_regularizer_weight 1e-4: Sets the regularization weight to 1e-4.
-tao_brgn_l1_smooth_epsilon 1e-6: Sets the smoothing parameter for L1 regularization to 1e-6.
-tao_gatol 1.e-5: Sets the gradient tolerance for convergence to 1.e-5.

4)   test:
      suffix: 4
      args: -tao_monitor_short -tao_max_it 100 -tao_type brgn -tao_brgn_regularization_type 
      lm -tao_gatol 1.e-5 -tao_brgn_subsolver_tao_type bnls

features:
-tao_monitor_short: Enables short monitoring of the TAO solver's progress.
-tao_max_it 100: Sets the maximum number of iterations to 100.
-tao_type brgn: Uses the BRGN (Bounded Regularized Gauss-Newton) solver type.
In this context, LM refers to the Levenberg-Marquardt algorithm, which is a popular 
method for solving nonlinear least squares problems. This version of BRGN combines 
bounded regularized Gauss-Newton with the damping strategy of Levenberg-Marquardt.
Key Features:
Adds a damping term to handle poorly conditioned problems or cases where the 
Gauss-Newton step might fail.
Balances between Gauss-Newton (for fast convergence) and gradient descent (for stability).
-tao_brgn_regularization_type lm: Uses Levenberg-Marquardt regularization.
-tao_gatol 1.e-5: Sets the gradient tolerance for convergence to 1.e-5.
-tao_brgn_subsolver_tao_type bnls: Uses the BNLS (Bounded Nonlinear Least Squares) 
                                    subsolver for the BRGN solver.

TEST*/