//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void SubStep1_y_cpu (real dt) {

//<USER_DEFINED>
  INPUT(Pressure);
  INPUT(Density);
  INPUT(Pot);
#ifdef X
  INPUT(Vx);
#ifdef COLLISIONPREDICTOR
  INPUT(Vx_half);
#endif
#endif
#ifdef Y
  INPUT(Vy);
#ifdef COLLISIONPREDICTOR
  INPUT(Vy_half);
#endif
  OUTPUT(Vy_temp);
#endif
#ifdef Z
#ifdef COLLISIONPREDICTOR
  INPUT(Vz_half);
#endif
  INPUT(Vz);
#endif
#ifdef MHD
  INPUT(Bx);
  INPUT(Bz);
#endif
  real xplanet = Xplanet;
  real yplanet = Yplanet; //para el sink
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
#ifdef COLLISIONPREDICTOR
  real* vx_half = Vx_half->field_cpu;
#else
  real* vx_half = Vx->field_cpu;
#endif
  real* vx_temp = Vx_temp->field_cpu;
#endif
#ifdef Y
  real* vy      = Vy->field_cpu;
#ifdef COLLISIONPREDICTOR
  real* vy_half = Vy_half->field_cpu;
#else
  real* vy_half = Vy->field_cpu;
#endif
  real* vy_temp = Vy_temp->field_cpu;
#endif
#ifdef Z
  real* vz      = Vz->field_cpu;
#ifdef COLLISIONPREDICTOR
  real* vz_half = Vz_half->field_cpu;
#else
  real* vz_half = Vz->field_cpu;
#endif
  real* vz_temp = Vz_temp->field_cpu;
#endif
#ifdef MHD
  real* bx = Bx->field_cpu;
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
  real delta = DELTA; // Parametros sink
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
#ifdef X
  int llxp;
#endif
#ifdef Y
  int llym;
#endif
#ifdef Z
  int llzp;
#endif
#ifdef MHD
  real db1;
  real db2;
  real bmeanm;
  real bmean;
#endif
#ifndef CARTESIAN
  real vphi;
#endif
#ifdef SHEARINGBOX
  real vm;
#endif
#ifdef SPHERICAL
  real vzz;
#endif
  real dx; 
  real dy;
  real dist2;
  real sink; // Debería llevar *??
  real planet_distance;
  real planet_angle;
  real vstary;
  real alpha;
  real sinkmomv;
  real r;
  real phi;
  real theta;
  real phi_p   = 0.0;
  real theta_p = M_PI/2.0;
  real rp      = 1.0;    
//<\INTERNAL>


//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real OMEGAFRAME(1);
// real OORTA(1);
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
#ifdef Y
	ll = l;
#ifdef X
	llxp = lxp;
#endif //ENDIF X
	llym = lym;
#ifdef Z
	llzp = lzp;
#endif //ENDIF Z
	dtOVERrhom = 2.0*dt/(rho[ll]*(ymin(j+1)-ymin(j))+rho[llym]*(ymin(j)-ymin(j-1)))*(ymed(j)-ymed(j-1));

	vy_temp[ll] = vy[ll];
	if(fluidtype != DUST) vy_temp[ll]-=  dtOVERrhom*(p[ll]-p[llym])/(ymed(j)-ymed(j-1));
	
#ifdef CARTESIAN
#ifdef SHEARINGBOX

	vm = 0.25*(vx_half[ll]+vx_half[llxp]+vx_half[llym]+vx_half[llxp-pitch]);	
	vy_temp[l] += dt*(2.0*OMEGAFRAME*vm + 2.0*SHEARPARAM*OMEGAFRAME*OMEGAFRAME*ymin(j));

#ifdef DRAGFORCE
	if (Fluidtype == GAS) {
	  vy_temp[ll] += 2*ASPECTRATIO*ASPECTRATIO*OMEGAFRAME*R0*dt;
	}
#endif
	
#endif //ENDIF SHEARINGBOX

#ifdef MHD
	if(fluidtype==GAS) {
#ifndef PASSIVEMHD
	  bmean  = 0.5*(bx[ll] + bx[llxp]);
	  bmeanm = 0.5*(bx[llym] + bx[llxp-pitch]);
	  
	  db1 = (bmean*bmean-bmeanm*bmeanm);
	  
	  bmean  = 0.5*(bz[ll] + bz[llzp]);
	  bmeanm = 0.5*(bz[llym] + bz[llym+stride]);
	  
	  db2 = (bmean*bmean-bmeanm*bmeanm); //grad(b^2/2)
	  
	  vy_temp[ll] -= .5*dtOVERrhom*(db1 + db2)/((ymed(j)-ymed(j-1))*MU0);
	  
#endif //ENDIF !PASSIVEMHD
	}
#endif //ENDIF MHD
#endif //ENDIF CARTESIAN
	  
#ifdef CYLINDRICAL
	vphi = .25*(vx_half[ll] + vx_half[llxp] + vx_half[llym] + vx_half[llxp-pitch]);
	vphi += ymin(j)*OMEGAFRAME;
	vy_temp[ll] += vphi*vphi/ymin(j)*dt;
	
#ifdef MHD
	if(fluidtype==GAS) {
#ifndef PASSIVEMHD
	bmean  = 0.5*(bx[ll] + bx[llxp]);
	bmeanm = 0.5*(bx[llym] + bx[llxp-pitch]);
	
	vy_temp[ll] -= dtOVERrhom*.25*(bmean+bmeanm)*(bmean+bmeanm)/(MU0*ymin(j));
	
	db1 = (bmean*bmean-bmeanm*bmeanm);
	
	bmean  = 0.5*(bz[ll] + bz[llzp]);
	bmeanm = 0.5*(bz[llym] + bz[llym+stride]);
	
	db2 = (bmean*bmean-bmeanm*bmeanm); //-grad((bphi^2+bz^2)/2)
	
	vy_temp[ll] -= .5*dtOVERrhom*(db1 + db2)/((ymed(j)-ymed(j-1))*MU0);
	
#endif //END !PASSIVEMHD
	}
#endif //END MHD
#endif //END CYLINDRICAL

#ifdef SPHERICAL
	vphi =  .25*(vx_half[ll] + vx_half[llxp] + vx_half[llym] + vx_half[llxp-pitch]);
	vphi += ymin(j)*sin(zmed(k))*OMEGAFRAME;
	vzz = .25*(vz_half[ll] + vz_half[llzp]  + vz_half[llym] + vz_half[llzp-pitch]);
	vy_temp[ll] += (vphi*vphi + vzz*vzz)/ymin(j)*dt;
#ifdef MHD
	if(fluidtype==GAS) {
#ifndef PASSIVEMHD
	  bmean  = 0.5*(bx[ll] + bx[llxp]);
	  bmeanm = 0.5*(bx[llym] + bx[llxp-pitch]);
	  
	  vy_temp[ll] -= dtOVERrhom*.25*(bmean+bmeanm)*(bmean+bmeanm)/(MU0*ymin(j));//Source term -B_phi^2/(mu_0\rho r)
	  
	  db1 = (bmean*bmean-bmeanm*bmeanm);
	  
	  bmean  = 0.5*(bz[ll] + bz[llzp]);
	  bmeanm = 0.5*(bz[llym] + bz[llym+stride]);
	  
	  vy_temp[ll] -= dtOVERrhom*.25*(bmean+bmeanm)*(bmean+bmeanm)/(MU0*ymin(j));//Source term -B_\theta^2/(mu_0\rho r)
	  
	  db2 = (bmean*bmean-bmeanm*bmeanm); //-grad((bphi^2+btheta^2)/2)
	  
	  vy_temp[ll] -= .5*dtOVERrhom*(db1 + db2)/((ymed(j)-ymed(j-1))*MU0);
	  
#endif //END !PASSIVEMHD
	}
#endif //END MHD
#endif //END SPHERICAL

#ifdef POTENTIAL
	  r     = ymin(j);
    phi   = xmed(i);
    theta = M_PI*0.5;
    vy_temp[ll] += -dt*G*( MSTAR/r/r + mp*taper*(r-rp*(cos(phi_p-phi)*sin(theta)*sin(theta_p)+cos(theta)*cos(theta_p)))/pow(-2.*r*rp*(cos(phi_p-phi)*sin(theta)*sin(theta_p)+cos(theta)*cos(theta_p))+rp*rp+r*r+smoothing,1.5));

#endif //ENDIF POTENTIAL

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
  alpha = planet_distance - xmed(i); 

  // Sink term
  if (planet_distance < rsink){
    sink = (1 - dist2/(rsink * rsink)) * (1 - dist2/(rsink * rsink));
  }
  else{
    sink = 0;
  }

 // starvector v*
  vstary = vy[ll]*(cos(alpha)*cos(alpha)+delta*sin(alpha)*sin(alpha)) + vphi*(1-delta)*(sin(alpha)*cos(alpha));


  sinkmomv = -gammasink * omegab * sink * vstary;
  
  vy_temp[ll] += sinkmomv * dt;
#endif //ENDIF SINKMOM
#endif //ENDIF Y


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
