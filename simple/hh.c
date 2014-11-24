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

static const FLOAT dt = 0.0025; // [msec]

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



FLOAT t;                   // [msec]
FLOAT hh_v;                // [mV]

FLOAT hh_n;
FLOAT hh_m;
FLOAT hh_h;


static FLOAT initialize()
{
  t = 0.0;                   // [msec]
  hh_v = 0.0;                // [mV]

  hh_n = 0.32;
  hh_m = 0.05;
  hh_h = 0.60;
}

int hh(FLOAT stoptime)
{
  const FLOAT hh_cm = 1.0;         // [muF/cm^2]

  const FLOAT hh_gk_max = 36.;    // [mS/cm^2]
  const FLOAT hh_gna_max = 120.;    // [mS/cm^2]
  const FLOAT hh_gm = 0.3;         // [mS/cm^3]

  const FLOAT hh_e_k = -12.0;      // [mV]
  const FLOAT hh_e_na = 115.0;     // [mV]
  const FLOAT hh_v_rest = 10.613;  // [mV]

  const FLOAT i_start =  50.;
  const FLOAT i_stop  = 175.;


  //FLOAT hh_n = calc_alpha_n(hh_v) / (calc_alpha_n(hh_v) + calc_beta_n(hh_v));
  //FLOAT hh_m = calc_alpha_m(hh_v) / (calc_alpha_m(hh_v) + calc_beta_m(hh_v));
  //FLOAT hh_h = calc_alpha_h(hh_v) / (calc_alpha_h(hh_v) + calc_beta_h(hh_v));


  for(t=0.0; t<stoptime; t+=dt)
    {
      FLOAT i_inj;
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

      a_n = calc_alpha_n(hh_v);
      a_m = calc_alpha_m(hh_v);
      a_h = calc_alpha_h(hh_v);
      b_n = calc_beta_n(hh_v);
      b_m = calc_beta_m(hh_v);
      b_h = calc_beta_h(hh_v);

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
      hh_n = hh_n + dt * (a_n * (1.0 - hh_n) - b_n * hh_n);
      hh_m = hh_m + dt * (a_m * (1.0 - hh_m) - b_m * hh_m);
      hh_h = hh_h + dt * (a_h * (1.0 - hh_h) - b_h * hh_h);

      i_k  = hh_gk_max  * hh_n * hh_n * hh_n * hh_n * (hh_e_k - hh_v);
      i_na = hh_gna_max * hh_m * hh_m * hh_m * hh_h * (hh_e_na - hh_v);
      i_m  = hh_gm * (hh_v_rest - hh_v);

      if(t > i_start && t < i_stop){
	i_inj = calc_i_inj(i_inj);
      }else{
	i_inj = 0.0;
      }
      
      hh_v = hh_v + dt / hh_cm * (i_k + i_na + i_m + i_inj);
      //hh_v = 0.0;

      //printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", t, a_n, b_n, a_m, b_m, a_h, b_h, hh_n, hh_m, hh_h, i_k, i_na, i_m, i_inj, hh_v);
      printf("%f %f %f\n", t, i_inj, hh_v);
    }

  return(0);
}

int main()
{
  //printf("Hodgkin-Huxley equation\n");

  hh(200);
  return(0);

}
