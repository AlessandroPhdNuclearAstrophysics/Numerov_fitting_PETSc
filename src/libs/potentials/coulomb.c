#include<gsl/gsl_sf_coulomb.h>
#include"coulomb.h"

void
coulomb_fg_(double *eta, double * x, double * lam_F, 
 double *fl, double *gl, double * flp, double * glp)
{

  gsl_sf_result F, Fp, G, Gp;
  double Fe, Ge;
  int k_G = 0;

  gsl_sf_coulomb_wave_FG_e(*eta, *x, *lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  
  *fl = F.val;
  *gl = G.val;

  *flp = Fp.val;
  *glp = Gp.val;
 
}
