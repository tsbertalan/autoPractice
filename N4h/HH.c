#include "auto_f2c.h"
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*   HHneurons :         HH-neurons with (chemical) synaptic coupling     */
/* --------------------------    by Sung Joon Moon (11/30/2006)    ------ */
/* --------------------------      moon@arnold.princeton.edu       ------ */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int func (integer ndim, const doublereal *y, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *f, doublereal *dfdu, doublereal *dfdp) {

  int		i, ind, N;
  doublereal	a_m, b_m, a_h, b_h, a_n, b_n,
		tau_m, tau_h, tau_n,
		m_inf, h_inf, n_inf,
		sum_s, C,
		V_Na, V_K, V_L, V_syn,	/* for now, excitatory case (V_syn>V_eq ~ -65mV) */
		g_Na, g_K, g_L,
		omega[ndim/5];		/* omega should be treated differently */
  
  /* Evaluates the algebraic equations or ODE right hand side */
  
  /* Input arguments : */
  /*      ndim   :   Dimension of the ODE system */
  /*      y      :   State variables */
  /*      icp    :   Array indicating the free parameter(s) */
  /*      par    :   Equation parameters; par[0] = I, par[1] = g_syn, par[2] = <tau>  */
  
  /* Values to be returned : */
  /*      f      :   ODE right hand side values; (V,m,h,n,s) for each neuron */
  
  /* Normally unused Jacobian arguments : IJAC, DFDU, DFDP (see manual) */
  
  C = 1.;
  V_Na = 50.; V_K = -77.; V_L = -54.5; V_syn = 30.;
  g_Na = 120.; g_K = 36.; g_L = .3;

  N = ndim/5;		/* number of neurons */

  for (i=0; i<N; i++) {	/* uniform distribution */
    omega[i] = -.1 + i*.2/(doublereal)(N-1);
  }

//printf("N = %d\n",N);
  sum_s = 0.;
  for (i=0; i<N; i++) {
    sum_s += y[i*5+4];
  }

  for (i=0; i<N; i++) {
    ind = i*5;		/* starting index for y-array, for the i-th neuron */

    a_m = .1*(y[ind]+40.)/(1.-exp((-y[ind]-40.)/10.));
    b_m = 4.*exp((-y[ind]-65.)/18.);
    a_h = .07*exp((-y[ind]-65.)/20.);
    b_h = 1./(1.+exp((-y[ind]-35.)/10.));
    a_n = .01*(y[ind]+55.)/(1.-exp((-y[ind]-55.)/10.));
    b_n = .125*exp((-y[ind]-65.)/80.);

    tau_m = 1./(a_m+b_m); m_inf = a_m*tau_m;
    tau_h = 1./(a_h+b_h); h_inf = a_h*tau_h;
    tau_n = 1./(a_n+b_n); n_inf = a_n*tau_n;

    f[ind]   = ( par[0] - g_Na*pow(y[ind+1],3.)*y[ind+2]*(y[ind]-V_Na)
		 - g_K*pow(y[ind+3],4.)*(y[ind]-V_K)
		 - g_L*(y[ind]-V_L)
		 - par[1]/N*(sum_s-y[ind+4])*(y[ind]-V_syn)
               )/C;
//printf("resid = %5.3e\n",f[ind]);
    f[ind+1] = (m_inf-y[ind+1])/tau_m;
    f[ind+2] = (h_inf-y[ind+2])/tau_h;
    f[ind+3] = (n_inf-y[ind+3])/tau_n;
    f[ind+4] = 1./(1.+exp(-(y[ind]-0.)/5.))*(1.-y[ind+4]) - y[ind+4]/(par[2]+omega[i]);
  }
  
  return 0;
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int stpnt (integer ndim, doublereal t,
           doublereal *y, doublereal *par) {
  /* Input arguments : */
  /*      ndim   :   Dimension of the ODE system */
  
  /* Values to be returned : */
  /*      y      :   A starting solution vector */
  /*      par    :   The corresponding equation-parameter values */
  
  
  /* Initialize the equation parameters */
  /*      par    :   Equation parameters; par[0] = I, par[1] = g_syn, par[2] = <tau>  */

  ndim = 20;

//  par[0] = (doublereal) 0.0;
  par[0] = (doublereal) 6.3931590417;
  par[1] = (doublereal) 3.0;
  par[2] = (doublereal) 1.0;
  
  /* Initialize the solution */
/*
  y[0] = -6.502507e+01;
  y[1] = 5.277622e-02;
  y[2] = 5.969972e-01;
  y[3] = 3.172928e-01;
  y[4] = 2.024113e-06;
  y[5] = -6.502508e+01;
  y[6] = 5.277616e-02;
  y[7] = 5.969976e-01;
  y[8] = 3.172926e-01;
  y[9] = 2.174043e-06;
  y[10] = -6.502509e+01;
  y[11] = 5.277611e-02;
  y[12] = 5.969979e-01;
  y[13] = 3.172925e-01;
  y[14] = 2.323973e-06;
  y[15] = -6.502510e+01;
  y[16] = 5.277605e-02;
  y[17] = 5.969982e-01;
  y[18] = 3.172924e-01;
  y[19] = 2.473902e-06;
*/
  return 0;
}
/* The following subroutines are not used here, */
/* but they must be supplied as dummy routines */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int bcnd (integer ndim, const doublereal *par, const integer *icp,
          integer nbc, const doublereal *u0, const doublereal *u1, integer ijac,
          doublereal *fb, doublereal *dbc) {
  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int icnd (integer ndim, const doublereal *par, const integer *icp,
          integer nint, const doublereal *u, const doublereal *uold,
          const doublereal *udot, const doublereal *upold, integer ijac,
          doublereal *fi, doublereal *dint) {
    return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int fopt (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *fs, doublereal *dfdu, doublereal *dfdp) {
    return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int pvls (integer ndim, const doublereal *u,
          doublereal *par) {
    return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
