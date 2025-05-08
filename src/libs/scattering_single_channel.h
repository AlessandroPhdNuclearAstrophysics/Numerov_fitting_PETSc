#ifndef SCATTERING_SINGLE_CHANNEL_ 
#define SCATTERING_SINGLE_CHANNEL_

#include<petsc.h>

/**
 * @struct coeffs
 * @brief Stores coefficients for scattering calculations.
 */
typedef struct {
  PetscReal cba[3]; /**< Array of three real numbers representing coefficients. */
} coeffs;

/**
 * @struct Observables
 * @brief Contains an array of coefficients for observables.
 */
typedef struct {
  coeffs* coeffs; /**< Pointer to an array of coeffs structures. */
} Observables;

/**
 * @struct QuantumNumbers
 * @brief Represents quantum numbers used in scattering calculations.
 */
typedef struct {
  PetscInt ipot; /**< Potential index. */
  PetscInt ilb;  /**< Lower bound index. */
  PetscInt L;    /**< Orbital angular momentum quantum number. */
  PetscInt S;    /**< Spin quantum number. */
  PetscInt *J;   /**< Pointer to an array of total angular momentum quantum numbers. */
  PetscInt nJ;   /**< Number of total angular momentum quantum numbers. */
} QuantumNumbers;

/**
 * @brief Computes scattering coefficients using the Numerov method.
 * 
 * @param energies Pointer to an array of energy values.
 * @param kcotd Pointer to an array for storing k*cot(delta) values.
 * @param ne Number of energy values.
 * @param ipot Potential index.
 * @param ilb Lower bound index.
 * @param L Orbital angular momentum quantum number.
 * @param S Spin quantum number.
 * @param J Total angular momentum quantum number.
 * @param r Radius parameter.
 * @param c Constant parameter.
 * @return A coeffs structure containing the computed coefficients.
 */
coeffs scattering_numerov(const double* energies, double *kcotd, const int ne,
        const int ipot, const int ilb, const int L, const int S, const int J, int nc, const double R, const double *LECS) {
  coeffs coeff;
  
  extern void numerov_scattering_single_channel_(
    const double *energies, 
    double *k_cot_delta, 
    const int *npoints_energy, 
    const int *L, 
    const int *S,
    const int *J, 
    const int *ipot, 
    const int *ilb, 
    const int *nc,
    const double *LECS, 
    const double *R, 
    double *coeff
  );

  if (LECS == NULL) {
    const double tmp[] = {0.0};
    numerov_scattering_single_channel_(energies, kcotd, &ne, &L, &S, &J, &ipot, &ilb, &nc, tmp, &R, (double *)&coeff);
    return coeff;
  }
  numerov_scattering_single_channel_(energies, kcotd, &ne, &L, &S, &J, &ipot, &ilb, &nc, LECS, &R, (double *)&coeff);
  return coeff;  
}

/**
 * @brief Computes observables for a given set of quantum numbers using the Numerov method.
 * 
 * @param energies Pointer to an array of energy values.
 * @param kcotd Array of pointers for storing k*cot(delta) values for each J.
 * @param ne Number of energy values.
 * @param qn Pointer to a QuantumNumbers structure containing quantum numbers.
 * @param r Radius parameter.
 * @param c Constant parameter.
 * @return An Observables structure containing the computed coefficients for each J.
 */
Observables scattering_numerov_quantum_num(const double* energies, double *kcotd[], const int ne,
                                                   const QuantumNumbers *qn, int nc, const double *LECS) {
  Observables observables;
  observables.coeffs = (coeffs*) malloc(qn->nJ * sizeof(coeffs));
  if (observables.coeffs == NULL) {
    fprintf(stderr, "Memory allocation failed for observables.coeffs\n");
    exit(EXIT_FAILURE);
  }
  int ipot = qn->ipot;
  int ilb = qn->ilb;
  int L = qn->L;
  int S = qn->S;
  
  for (int i = 0; i < qn->nJ; i++) {
    int J = qn->J[i];
    observables.coeffs[i] = scattering_numerov(energies, kcotd[i], ne, ipot, ilb, L, S, J, nc, LECS[0], LECS);
  }
  return observables;
}

#endif // SCATTERING_SINGLE_CHANNEL_