#include <stdio.h>
#include <math.h>

typedef double FLOAT;
// typedef double FLOAT;
#define EXP(x) exp(x)
// #define EXP(x) expf(x)

FLOAT alpha_m (FLOAT v) { return( (25. - v)/(10. * EXP( (25. - v) / 10.) - 1.) ); }
FLOAT beta_m  (FLOAT v) { return( 4. * EXP(- v / 18.));                           }
FLOAT alpha_h (FLOAT v) { return( 0.07 * EXP(- v / 20.) );                        }
FLOAT beta_h  (FLOAT v) { return( 1. / (EXP(30. - v) / 10. + 1.) );               }
FLOAT alpha_n (FLOAT v) { return( (10. - v)/(100. * EXP((10. - v) / 10.) - 1.) ); }
FLOAT beta_n  (FLOAT v) { return( 0.125 * EXP(- v / 80.));                        }

int hh(FLOAT stoptime)
{
  FLOAT dt = 0.025;          // [msec]

  FLOAT hh_cm = 1.0;         // [muF/cm^2]

  FLOAT hh_gna_max = 120.0;  // [mS/cm^2]
  FLOAT hh_gk_max = 36.0;    // [mS/cm^2]
  FLOAT hh_gm = 0.3;         // [mS/cm^3]

  FLOAT hh_e_na = 115.0;     // [mV]
  FLOAT hh_e_k = -12.0;      // [mV]

  FLOAT hh_v_rest = 10.613;  // [mV]
  FLOAT hh_v = 0.0;          // [mV]

  FLOAT hh_m = 0.05;
  FLOAT hh_h = 0.60;
  FLOAT hh_n = 0.32;


  FLOAT i_na  = 0.0;
  FLOAT i_k   = 0.0;
  FLOAT i_m   = 0.0;
  FLOAT i_inj = 0.0;


}

int main()
{
  printf("Hodgkin-Huxley equation\n");

  return(0);

}
