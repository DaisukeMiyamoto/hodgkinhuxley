#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef double FLOAT;
// typedef float FLOAT;
//#define EXP(x) exp((x))
// #define EXP(x) expf(x)
#define EXP exp

static FLOAT calc_alpha_n (FLOAT v) { return( (10. - v) / (100. * (EXP( (10. - v) / 10.) - 1.) ) ); }
static FLOAT calc_beta_n  (FLOAT v) { return( 0.125 * EXP(- v / 80.));                         }
static FLOAT calc_alpha_m (FLOAT v) { return( (25. - v) / (10. * (EXP( (25. - v) / 10.) - 1.) ) );  }
static FLOAT calc_beta_m  (FLOAT v) { return( 4. * EXP(- v / 18.) );                            }
static FLOAT calc_alpha_h (FLOAT v) { return( 0.07 * EXP(- v / 20.) );                         }
static FLOAT calc_beta_h  (FLOAT v) { return( 1. / (EXP( (30. - v) / 10.) + 1.) );             }

static const FLOAT dt = 0.025; // [msec]

static FLOAT calc_i_inj(FLOAT i_inj)
{
  const FLOAT i_dc   = 10.;
  const FLOAT i_rand = 0.;
  const FLOAT i_tau  = 1.;

  FLOAT i_inj_raw;

  i_inj_raw = i_dc + ((float) rand() / RAND_MAX - 0.5) * 2. * i_rand;
  i_inj  += dt / i_tau * (i_inj_raw - i_inj);

  return(i_inj);
}

#define N_COMPARTMENT 10000

FLOAT hh_v[N_COMPARTMENT];                // [mV]
FLOAT hh_n[N_COMPARTMENT];
FLOAT hh_m[N_COMPARTMENT];
FLOAT hh_h[N_COMPARTMENT];

const FLOAT hh_cm[N_COMPARTMENT]      = {1.0};    // [muF/cm^2]
const FLOAT hh_gk_max[N_COMPARTMENT]  = {36.};    // [mS/cm^2]
const FLOAT hh_gna_max[N_COMPARTMENT] = {120.};   // [mS/cm^2]
const FLOAT hh_gm[N_COMPARTMENT]      = {0.3};    // [mS/cm^3]

const FLOAT hh_e_k        = -12.0;    // [mV]
const FLOAT hh_e_na       = 115.0;    // [mV]
const FLOAT hh_v_rest     = 10.613;   // [mV]


static void initialize()
{
  int i;
  for(i=0; i<N_COMPARTMENT; i++)
    {
      hh_v[i] = 0.0;                // [mV]
      hh_n[i] = 0.32;
      hh_m[i] = 0.05;
      hh_h[i] = 0.60;
      //FLOAT hh_n = calc_alpha_n(hh_v) / (calc_alpha_n(hh_v) + calc_beta_n(hh_v));
      //FLOAT hh_m = calc_alpha_m(hh_v) / (calc_alpha_m(hh_v) + calc_beta_m(hh_v));
      //FLOAT hh_h = calc_alpha_h(hh_v) / (calc_alpha_h(hh_v) + calc_beta_h(hh_v));
    }
  return;
}

int hh(FLOAT stoptime)
{
  int i, j;
  int i_stop;

  const int inj_start =  50./dt;
  const int inj_stop  = 175./dt;

  initialize();
  for(i=0,i_stop=stoptime/dt; i<i_stop; i++)
    {
      FLOAT i_inj;
      if(i > inj_start && i < inj_stop){
	i_inj = calc_i_inj(i_inj);
      }else{
	i_inj = 0.0;
      }

      for(j=0; j<N_COMPARTMENT; j++)
	{
	  FLOAT i_k;
	  FLOAT i_na;
	  FLOAT i_m;
	  
	  FLOAT a_n;
	  FLOAT a_m;
	  FLOAT a_h;
	  FLOAT b_n;
	  FLOAT b_m;
	  FLOAT b_h;
	  FLOAT tau_n;
	  FLOAT tau_m;
	  FLOAT tau_h;
	  FLOAT n_inf;
	  FLOAT m_inf;
	  FLOAT h_inf;
	  
	  a_n = calc_alpha_n(hh_v[j]);
	  a_m = calc_alpha_m(hh_v[j]);
	  a_h = calc_alpha_h(hh_v[j]);
	  b_n = calc_beta_n(hh_v[j]);
	  b_m = calc_beta_m(hh_v[j]);
	  b_h = calc_beta_h(hh_v[j]);
	  
	  tau_n = 1. / (a_n + b_n);
	  tau_m = 1. / (a_m + b_m);
	  tau_h = 1. / (a_h + b_h);
	  n_inf = a_n * tau_n;
	  m_inf = a_m * tau_m;
	  h_inf = a_h * tau_h;
	  
	  /*
	    hh_n += dt * (n_inf - hh_n) / tau_n;
	    hh_m += dt * (m_inf - hh_m) / tau_m;
	    hh_h += dt * (h_inf - hh_h) / tau_h;
	  */
	  /*
	    hh_n = hh_n + dt * (a_n * (1.0 - hh_n) - b_n * hh_n);
	    hh_m = hh_m + dt * (a_m * (1.0 - hh_m) - b_m * hh_m);
	    hh_h = hh_h + dt * (a_h * (1.0 - hh_h) - b_h * hh_h);
	  */
	  hh_n[j] += (1.0 - expf(-dt / tau_n)) * (n_inf - hh_n[j]);
	  hh_m[j] += (1.0 - expf(-dt / tau_m)) * (m_inf - hh_m[j]);
	  hh_h[j] += (1.0 - expf(-dt / tau_h)) * (h_inf - hh_h[j]);
	  
	  
	  i_k  = hh_gk_max[0]  * hh_n[j] * hh_n[j] * hh_n[j] * hh_n[j] * (hh_e_k - hh_v[j]);
	  i_na = hh_gna_max[0] * hh_m[j] * hh_m[j] * hh_m[j] * hh_h[j] * (hh_e_na - hh_v[j]);
	  i_m  = hh_gm[0] * (hh_v_rest - hh_v[j]);
	  
	  hh_v[j] = hh_v[j] + dt / hh_cm[0] * (i_k + i_na + i_m + i_inj);  
	}
      //printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", t, a_n, b_n, a_m, b_m, a_h, b_h, hh_n, hh_m, hh_h, i_k, i_na, i_m, i_inj, hh_v);
      printf("%f %f %f %f\n", i*dt, i_inj, hh_v[0], hh_v[20]);
    }

  return(0);
}

int main()
{
  //printf("Hodgkin-Huxley equation\n");

  hh(200);
  return(0);

}
