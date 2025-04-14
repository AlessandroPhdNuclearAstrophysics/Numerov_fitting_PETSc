#ifndef SCATTERING_SINGLE_CHANNEL_ 
#define SCATTERING_SINGLE_CHANNEL_

#include<petsc.h>

typedef struct {
  PetscReal cba[3];
} coeffs;

coeffs scattering_numerov(const double* energies, double *kcotd, const int ne, 
                          const int ipot, const int ilb, const double r, const double c) {
  coeffs coeff;
  int L=1, S=0, J=1;

  extern void numerov_scattering_single_channel_(const double *, double *, const int *, const int *, const int *,
    const int *, const int *, const int *, const double *, const double *, double *);

  numerov_scattering_single_channel_(energies, kcotd, &ne, &L, &S, &J, &ipot, &ilb, &c, &r, (double *)&coeff);
  return coeff;  
}


#endif