//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void SubStep1_x_cpu (real dt) {


//<USER_DEFINED>
  INPUT(Pressure);
  INPUT(Density);
  INPUT(Pot);
  INPUT(Vx);
  INPUT(Vy);
#ifdef MHD
#if defined (CYLINDRICAL) || defined (SPHERICAL)
  INPUT(Bx);
#endif
  INPUT(By);
  INPUT(Bz);
#endif  
  OUTPUT(Vx_temp);
  real xplanet = Xplanet;
  real yplanet = Yplanet;
  // potencial modif
  real planetmass_taper;
  if (MASSTAPER == 0.0)
    planetmass_taper = 1.0;
  else
    planetmass_taper = (PhysicalTime >= MASSTAPER ? 1.0 : .5*(1.0-cos(M_PI*PhysicalTime/MASSTAPER)));
//<\USER_DEFINED>

//<EXTERNAL>
  real* p   = Pressure->field_cpu;
  real* pot = Pot->field_cpu;
  real* rho = Density->field_cpu;
#ifdef X
  real* vx      = Vx->field_cpu;
  real* vx_temp = Vx_temp->field_cpu;
#endif
#ifdef Y 
  real* vy = Vy->field_cpu;
#endif
#ifdef MHD
  real* bx = Bx->field_cpu;
  real* by = By->field_cpu;
  real* bz = Bz->field_cpu;
#endif
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = XIP; 
  int size_y = Ny+2*NGHY-1;
  int size_z = Nz+2*NGHZ-1;
  int fluidtype = Fluidtype;
  real gammasink = GAMMASINK;
  real rsink = RSINK;
  real omegab = OMEGAB;
  real delta = DELTA;
  // potencial modif
  real mp = PLANETMASS;
  real smoothing = THICKNESSSMOOTHING*THICKNESSSMOOTHING;
  real taper = planetmass_taper;
//<\EXTERNAL>
  
//<INTERNAL>
  int i; //Variables reserved
  int j; //for the topology
  int k; //of the kernels
  int ll;
  real dtOVERrhom;
  real dxmed;
#ifdef X
  int llxm;
#endif
#ifdef Y
  int llyp;
#endif
#ifdef Z
  int llzp;
#endif
#ifdef MHD
  real db1;
  real db2;
  real bmeanm;
  real bmean;
#if defined (CYLINDRICAL)|| defined (SPHERICAL)
  real brmean;
#endif
#ifdef SPHERICAL
  real btmean;
#endif
#endif
  real dx; 
  real dy;
  real dist2;
  real sink;
  real sinkmomv;
  real planet_distance;
  real planet_angle;
  real vstarx;
  real alpha;
  real r;
  real phi;
  real theta;
  real phi_p   = 0.0;
  real theta_p = M_PI/2.0;
  real rp      = 1.0; 
 //<\INTERNAL>

//<CONSTANT>
// real Xplanet(1);
// real Yplanet(1);
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real InvDiffXmed(Nx+2*NGHX);
// real Sxi(Nx+2*NGHX);
//<\CONSTANT>


//<MAIN_LOOP>

  i = j = k = 0;
  
#ifdef Z
  for(k=1; k<size_z; k++) {
#endif
#ifdef Y
    for(j=1; j<size_y; j++) {
#endif
#ifdef X
      for(i=XIM; i<size_x; i++) {
#endif
//<#>
#ifdef X
	ll = l;
	llxm = lxm;
#ifdef Y
	llyp = lyp;
#endif
#ifdef Z
	llzp = lzp;
#endif

	dxmed = 0.5*( Sxi(i) + Sxi(ixm) );
	dtOVERrhom = 2.*dt/( rho[ll]*Sxi(i) + rho[llxm]*Sxi(ixm) )*dxmed;

	vx_temp[ll] = vx[ll];
	if(fluidtype != DUST) vx_temp[ll] -=  dtOVERrhom*(p[ll]-p[llxm])*Inv_zone_size_xmed(i,j,k);
	
#ifdef POTENTIAL
  r = ymed(j);
	phi = xmin(i);
	theta = M_PI*0.5;

	vx_temp[ll] += dt/(r*sin(theta))*G*mp*taper*((sin(phi_p-phi)*r*rp*sin(theta)*sin(theta_p))/pow(-2.*r*rp*(cos(phi_p-phi)*sin(theta)*sin(theta_p)+cos(theta)*cos(theta_p))+rp*rp+r*r+smoothing,1.5));

#endif

#ifdef SINKMOM
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
  planet_angle = atan(dy/dx); 
  alpha = planet_angle - xmed(i); 

  // Sink term
  if (planet_distance < rsink){
    sink = (1 - dist2/(rsink * rsink)) * (1 - dist2/(rsink * rsink));
  }
  else{
    sink = 0;
  }
  //
  // starvector v* 
 
  vstarx = 0.25*(vy[ll] + vy[lxm] + vy[lxm-pitch] + vy[lyp])*(1-delta)*(sin(alpha)*cos(alpha)) + vx[ll]*(sin(alpha)*sin(alpha)+delta*cos(alpha)*cos(alpha));
  
  //

  sinkmomv = -gammasink * omegab * sink * vstarx;
  
  vx_temp[ll] += sinkmomv * dt;
#endif

#ifdef MHD
	if(fluidtype == GAS) {

#ifndef PASSIVEMHD
	  
	  bmean  = 0.5*(by[ll] + by[llyp]);
	  bmeanm = 0.5*(by[llxm] + by[llxm+pitch]);
	  
	  db1 = (bmean*bmean-bmeanm*bmeanm);
	  
#if defined(CYLINDRICAL) || defined(SPHERICAL)
	  brmean = .5*(bmean+bmeanm);
	  vx_temp[ll] += dtOVERrhom*brmean*bx[ll]/(MU0*ymed(j));
#endif
	  
	  bmean  = 0.5*(bz[ll] + bz[llzp]);
	  bmeanm = 0.5*(bz[llxm] + bz[llxm+stride]);
	  
	  db2 = (bmean*bmean-bmeanm*bmeanm);
	  
	  vx_temp[ll] -= .5*dtOVERrhom*(db1 + db2)*Inv_zone_size_xmed(i,j,k)/MU0;
	  
#ifdef SPHERICAL
	btmean = .5*(bmean+bmeanm);
	vx_temp[ll] += dtOVERrhom*btmean*cos(zmed(k))*bx[ll]/(MU0*ymed(j)*sin(zmed(k)));
#endif
      
#endif
      }
#endif
#endif
      

//<\#>
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
//<\MAIN_LOOP>





}








