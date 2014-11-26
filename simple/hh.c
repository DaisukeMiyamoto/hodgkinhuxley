#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//#define EXP(x) exp((x))
// #define EXP(x) expf(x)
#define EXP exp

#define N_COMPARTMENT 10000

typedef double FLOAT;
// typedef float FLOAT;

static FLOAT calc_alpha_n (FLOAT v) { return( (10. - v) / (100. * (EXP( (10. - v) / 10.) - 1.) ) ); }
static FLOAT calc_beta_n  (FLOAT v) { return( 0.125 * EXP(- v / 80.));                         }
static FLOAT calc_alpha_m (FLOAT v) { return( (25. - v) / (10. * (EXP( (25. - v) / 10.) - 1.) ) );  }
static FLOAT calc_beta_m  (FLOAT v) { return( 4. * EXP(- v / 18.) );                            }
static FLOAT calc_alpha_h (FLOAT v) { return( 0.07 * EXP(- v / 20.) );                         }
static FLOAT calc_beta_h  (FLOAT v) { return( 1. / (EXP( (30. - v) / 10.) + 1.) );             }


static const FLOAT dt = 0.025; // [msec]

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


#define TABLE_N_TAU(x) hh_table[(x)][0]
#define TABLE_N_INF(x) hh_table[(x)][1]
#define TABLE_M_TAU(x) hh_table[(x)][2]
#define TABLE_M_INF(x) hh_table[(x)][3]
#define TABLE_H_TAU(x) hh_table[(x)][4]
#define TABLE_H_INF(x) hh_table[(x)][5]
#define TABLE_SIZE 201
#define TABLE_MAX_V 150.0
#define TABLE_MIN_V -50.0

FLOAT hh_table[TABLE_SIZE][6];

static void makeTable()
{
  int i;
  for(i=0; i<TABLE_SIZE; i++)
    {
      FLOAT v;
      FLOAT a_n, a_m, a_h, b_n, b_m, b_h;
      v = (TABLE_MAX_V - TABLE_MIN_V)/(FLOAT)TABLE_SIZE * i + TABLE_MIN_V;
      a_n = calc_alpha_n(v);
      a_m = calc_alpha_m(v);
      a_h = calc_alpha_h(v);
      b_n = calc_beta_n(v);
      b_m = calc_beta_m(v);
      b_h = calc_beta_h(v);

      TABLE_N_TAU(i) = 1. / (a_n + b_n);
      TABLE_N_INF(i) = a_n * TABLE_N_TAU(i);
      TABLE_M_TAU(i) = 1. / (a_m + b_m);
      TABLE_M_INF(i) = a_m * TABLE_M_TAU(i);
      TABLE_H_TAU(i) = 1. / (a_h + b_h);
      TABLE_H_INF(i) = a_h * TABLE_H_TAU(i);
    }
  return;
}


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
	  hh_n[j] += (1.0 - EXP(-dt / tau_n)) * (n_inf - hh_n[j]);
	  hh_m[j] += (1.0 - EXP(-dt / tau_m)) * (m_inf - hh_m[j]);
	  hh_h[j] += (1.0 - EXP(-dt / tau_h)) * (h_inf - hh_h[j]);
	  
	  
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

int hh_with_table(FLOAT stoptime)
{
  int i, j;
  int i_stop;
  const int inj_start =  50./dt;
  const int inj_stop  = 175./dt;

  initialize();
  makeTable();

  for(i=0,i_stop=stoptime/dt; i<i_stop; i++)
    {
      FLOAT i_inj;
      if(i > inj_start && i < inj_stop){
	i_inj = 10.0;
      }else{
	i_inj = 0.0;
      }

      for(j=0; j<N_COMPARTMENT; j++)
	{
	  FLOAT i_k, i_na, i_m;
	  FLOAT tau_n, n_inf, tau_m, m_inf, tau_h, h_inf;
	  FLOAT theta;
	  int v_i;
	  
	  v_i = (int)(hh_v[j] - TABLE_MIN_V);
	  theta = (hh_v[j] - TABLE_MIN_V) - (FLOAT)v_i;
	  if(v_i>=TABLE_SIZE){ v_i=TABLE_SIZE-1; theta=1.0; }
	  if(v_i<0){ v_i=0; theta=0.0; }

	  tau_n = TABLE_N_TAU(v_i) + theta * (TABLE_N_TAU(v_i+1) - TABLE_N_TAU(v_i));
	  n_inf = TABLE_N_INF(v_i) + theta * (TABLE_N_INF(v_i+1) - TABLE_N_INF(v_i));
	  tau_m = TABLE_M_TAU(v_i) + theta * (TABLE_M_TAU(v_i+1) - TABLE_M_TAU(v_i));
	  m_inf = TABLE_M_INF(v_i) + theta * (TABLE_M_INF(v_i+1) - TABLE_M_INF(v_i));
	  tau_h = TABLE_H_TAU(v_i) + theta * (TABLE_H_TAU(v_i+1) - TABLE_H_TAU(v_i));
	  h_inf = TABLE_H_INF(v_i) + theta * (TABLE_H_INF(v_i+1) - TABLE_H_INF(v_i));

	  hh_n[j] += (1.0 - EXP(-dt / tau_n)) * (n_inf - hh_n[j]);
	  hh_m[j] += (1.0 - EXP(-dt / tau_m)) * (m_inf - hh_m[j]);
	  hh_h[j] += (1.0 - EXP(-dt / tau_h)) * (h_inf - hh_h[j]);	  
	  
	  i_k  = hh_gk_max[0]  * hh_n[j] * hh_n[j] * hh_n[j] * hh_n[j] * (hh_e_k - hh_v[j]);
	  i_na = hh_gna_max[0] * hh_m[j] * hh_m[j] * hh_m[j] * hh_h[j] * (hh_e_na - hh_v[j]);
	  i_m  = hh_gm[0] * (hh_v_rest - hh_v[j]);
	  
	  hh_v[j] = hh_v[j] + dt / hh_cm[0] * (i_k + i_na + i_m + i_inj);  
	  
	}
      printf("%f %f %f %f\n", i*dt, i_inj, hh_v[0], hh_v[20]);
    }

  return(0);
}



int main()
{
  //printf("Hodgkin-Huxley equation\n");

  //hh(200);
  hh_with_table(200);
  return(0);

}
