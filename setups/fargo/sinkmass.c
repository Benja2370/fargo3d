//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void Sinkdens(real dt) {

//<USER_DEFINED>
  INPUT(Density);
  OUTPUT(Density);
  // Para planeta fijo desde potencial.c
  real xplanet = 1.0;
  real yplanet = 0.0; //para el sink
//<\USER_DEFINED>

//<EXTERNAL>
  real* dens = Density->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  real gammasink = GAMMASINK;
  real rsink = RSINK;
  real omegab = OMEGAB;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  real dx;
  real dy;
  real dist2;
  real planet_distance;
  real sink;
  real sinkdens;
//<\INTERNAL>


//<MAIN_LOOP>
  for (k=0; k<size_z; k++) {
    for (j=0; j<size_y; j++) {
      for (i=0; i<size_x; i++) {
//<#>   
        #ifdef CYLINDRICAL
	        dx = ymed(j)*cos(xmed(i))-xplanet;
	        dy = ymed(j)*sin(xmed(i))-yplanet;
        #endif
        #ifdef SPHERICAL
	        dx = ymed(j)*cos(xmed(i))*sin(zmed(k))-xplanet;
	        dy = ymed(j)*sin(xmed(i))*sin(zmed(k))-yplanet;
        #endif

        dist2 = dx*dx+dy*dy;
        planet_distance = sqrt(dist2);
        
        // Termino sink
        if (planet_distance < rsink){
            sink = (1 - dist2/(rsink * rsink)) * (1 - dist2/(rsink * rsink));
        }
        else{
            sink = 0;
        }
        //
        sinkdens = -gammasink * omegab * sink * dens[l]; //interpolación?

        dens[l] += sinkdens * dt;

//<\#>
      }
    }
  }
//<\MAIN_LOOP>
}
