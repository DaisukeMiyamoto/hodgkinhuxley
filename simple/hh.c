#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//#define EXP(x) exp((x))
// #define EXP(x) expf(x)
#define EXP exp

#define N_COMPARTMENT 10000

typedef double FLOAT;
//typedef float FLOAT;


static const FLOAT DT = 0.025;        // [msec]
static const FLOAT CELSIUS = 6.3;     // []

/*
static FLOAT calc_alpha_n (FLOAT v) { return( (10. - v) / (100. * (EXP( (10. - v) / 10.) - 1.) ) ); }
static FLOAT calc_beta_n  (FLOAT v) { return( 0.125 * EXP(- v / 80.));                         }
static FLOAT calc_alpha_m (FLOAT v) { return( (25. - v) / (10. * (EXP( (25. - v) / 10.) - 1.) ) );  }
static FLOAT calc_beta_m  (FLOAT v) { return( 4. * EXP(- v / 18.) );                            }
static FLOAT calc_alpha_h (FLOAT v) { return( 0.07 * EXP(- v / 20.) );                         }
static FLOAT calc_beta_h  (FLOAT v) { return( 1. / (EXP( (30. - v) / 10.) + 1.) );             }
*/
static FLOAT calc_alpha_n (FLOAT v) { return( 0.01  * -(v+55.) / (EXP( -(v+55.)/10.) - 1.) ); }
static FLOAT calc_beta_n  (FLOAT v) { return( 0.125 * EXP( -(v+65.) / 80.) );                 }
static FLOAT calc_alpha_m (FLOAT v) { return( 0.1   * -(v+40.) / (EXP( -(v+40.)/10.) - 1.) ); }
static FLOAT calc_beta_m  (FLOAT v) { return( 4.    * EXP( -(v+65) / 18.) );                  }
static FLOAT calc_alpha_h (FLOAT v) { return( 0.07  * EXP( -(v+65) / 20.) );                  }
static FLOAT calc_beta_h  (FLOAT v) { return( 1. / (EXP( -(v+35) / 10.) + 1.) );              }

//static FLOAT correct_temp (void) { return(pow(3.0, CELSIUS - 6.3 / 10.)); }
//static FLOAT calc_n_inf (FLOAT v) { return( 1.0 / (correct_temp() * (calc_alpha_n(v)+calc_beta_n(v)) ) ); }


static FLOAT hh_v[N_COMPARTMENT];     // [mV]
static FLOAT hh_n[N_COMPARTMENT];
static FLOAT hh_m[N_COMPARTMENT];
static FLOAT hh_h[N_COMPARTMENT];


static FLOAT hh_cm[N_COMPARTMENT];
static FLOAT hh_cm_inv[N_COMPARTMENT];
static FLOAT hh_gk_max[N_COMPARTMENT];
static FLOAT hh_gna_max[N_COMPARTMENT];
static FLOAT hh_gm[N_COMPARTMENT];

const FLOAT hh_e_k        = -77.0;    // [mV]
const FLOAT hh_e_na       =  50.0;    // [mV]
const FLOAT hh_v_rest     = -54.3;    // [mV]


static void initialize()
{
  int i;
  for(i=0; i<N_COMPARTMENT; i++)
    {
      hh_v[i] = -65.0;                // [mV]

      //hh_n[i] = 0.32;
      //hh_m[i] = 0.05;
      //hh_h[i] = 0.60;
      hh_n[i] = calc_alpha_n(hh_v[i]) / (calc_alpha_n(hh_v[i]) + calc_beta_n(hh_v[i]));
      hh_m[i] = calc_alpha_m(hh_v[i]) / (calc_alpha_m(hh_v[i]) + calc_beta_m(hh_v[i]));
      hh_h[i] = calc_alpha_h(hh_v[i]) / (calc_alpha_h(hh_v[i]) + calc_beta_h(hh_v[i]));

      hh_cm[i]      = 1.0;             // [muF/cm^2]
      hh_cm_inv[i]  = 1.0 / hh_cm[i];  // [cm^2/muF]
      hh_gk_max[i]  = 36.;           // [S/cm^2]
      hh_gna_max[i] = 120.;            // [S/cm^2] 
      hh_gm[i]      = 0.3;          // [S/cm^3]
    }

  /*
  FLOAT v=-65.;
  FLOAT a_n, a_m, a_h, b_n, b_m, b_h;
  a_n = calc_alpha_n(v);
  a_m = calc_alpha_m(v);
  a_h = calc_alpha_h(v);
  b_n = calc_beta_n(v);
  b_m = calc_beta_m(v);
  b_h = calc_beta_h(v);
  printf("n (tau, inf)=(%.4f, %.4f)\n", 1./(a_n + b_n), a_n / (a_n+b_n));
  printf("m (tau, inf)=(%.4f, %.4f)\n", 1./(a_m + b_m), a_m / (a_m+b_m));
  printf("h (tau, inf)=(%.4f, %.4f)\n", 1./(a_h + b_h), a_h / (a_h+b_h));
  */

  return;
}


#define TABLE_N_TAU(x) hh_table[(x)][0]
#define TABLE_N_INF(x) hh_table[(x)][1]
#define TABLE_M_TAU(x) hh_table[(x)][2]
#define TABLE_M_INF(x) hh_table[(x)][3]
#define TABLE_H_TAU(x) hh_table[(x)][4]
#define TABLE_H_INF(x) hh_table[(x)][5]
#define TABLE_SIZE 201
#define TABLE_MAX_V 100.0
#define TABLE_MIN_V -100.0

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
  i_inj  += DT / i_tau * (i_inj_raw - i_inj);

  return(i_inj);
}


int hh_with_table(FLOAT stoptime)
{
  unsigned int i, j;
  unsigned int i_stop;
  unsigned int v_i_array[N_COMPARTMENT];
  FLOAT theta_array[N_COMPARTMENT];

  const int inj_start =  50./DT;
  const int inj_stop  = 175./DT;

  initialize();
  makeTable();

  for(i=0,i_stop=stoptime/DT; i<i_stop; i++)
    {
      FLOAT i_inj;
      if(i > inj_start && i < inj_stop){
	i_inj = 10.0 / 1000.;
      }else{
	i_inj = 0.0;
      }

      for(j=0; j<N_COMPARTMENT; j++)
	{
	  v_i_array[j] = (int)(hh_v[j] - TABLE_MIN_V);
	  theta_array[j] = (hh_v[j] - TABLE_MIN_V) - (FLOAT)v_i_array[j];
	  if(v_i_array[j] >= TABLE_SIZE){ v_i_array[j]=TABLE_SIZE-1; theta_array[j]=1.0; }
	  if(v_i_array[j] <  0)         { v_i_array[j]=0;            theta_array[j]=0.0; }
	}

      for(j=0; j<N_COMPARTMENT; j++)
	{
	  FLOAT tau_n, n_inf, tau_m, m_inf, tau_h, h_inf;
	  unsigned int v_i = v_i_array[j];
	  FLOAT theta = theta_array[j];

	  tau_n = TABLE_N_TAU(v_i) + theta * (TABLE_N_TAU(v_i+1) - TABLE_N_TAU(v_i));
	  n_inf = TABLE_N_INF(v_i) + theta * (TABLE_N_INF(v_i+1) - TABLE_N_INF(v_i));
	  tau_m = TABLE_M_TAU(v_i) + theta * (TABLE_M_TAU(v_i+1) - TABLE_M_TAU(v_i));
	  m_inf = TABLE_M_INF(v_i) + theta * (TABLE_M_INF(v_i+1) - TABLE_M_INF(v_i));
	  tau_h = TABLE_H_TAU(v_i) + theta * (TABLE_H_TAU(v_i+1) - TABLE_H_TAU(v_i));
	  h_inf = TABLE_H_INF(v_i) + theta * (TABLE_H_INF(v_i+1) - TABLE_H_INF(v_i));

	  hh_n[j] += (1.0 - EXP(-DT / tau_n)) * (n_inf - hh_n[j]);
	  hh_m[j] += (1.0 - EXP(-DT / tau_m)) * (m_inf - hh_m[j]);
	  hh_h[j] += (1.0 - EXP(-DT / tau_h)) * (h_inf - hh_h[j]);	  
	}
	  
      for(j=0; j<N_COMPARTMENT; j++)
	{
	  FLOAT i_k, i_na, i_m;
	  i_k  = hh_gk_max[i]  * hh_n[j] * hh_n[j] * hh_n[j] * hh_n[j] * (hh_e_k - hh_v[j]);
	  i_na = hh_gna_max[i] * hh_m[j] * hh_m[j] * hh_m[j] * hh_h[j] * (hh_e_na - hh_v[j]);
	  i_m  = hh_gm[i] * (hh_v_rest - hh_v[j]);
	  
	  hh_v[j] += DT * hh_cm_inv[i] * (i_k + i_na + i_m + i_inj);  
	  
	}
      printf("%f %f %f %f\n", i*dt, i_inj, hh_v[0], hh_v[N_COMPARTMENT-1]);
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
