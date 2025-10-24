/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Disk-Planet interaction problem.

*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include <stdbool.h>

#define MIN_DENSITY 1e-10
#if ROTATING_FRAME == NO
 #define g_OmegaZ  0.0
#endif

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double r, th, R, z, H, OmegaK, Omega, cs, cs02;
  double sigma;
  double p_index, q_index;
  double unit_surface_density = CONST_Msun/(UNIT_LENGTH*UNIT_LENGTH);

  p_index = -g_inputParam[P_INDEX];
  q_index = -g_inputParam[Q_INDEX];

 #if GEOMETRY == POLAR
  R = x1;
 #elif GEOMETRY == SPHERICAL
  R = x1*sin(x2);
 #endif

  sigma = (g_inputParam[SIGMA_0] / unit_surface_density) * pow(R, p_index);

  OmegaK = sqrt(G_MU/R)/R;

  // Calculate sound speed
#if FIXED_ASPECT_RATIO == NO
  double temperature = g_inputParam[TEMPERATURE_0] * pow(R, q_index);
  cs = sqrt(temperature/KELVIN/g_inputParam[MEAN_MOLECULAR_WEIGTH]);
 #if EOS == ISOTHERMAL
  cs02 = g_inputParam[TEMPERATURE_0]/KELVIN/g_inputParam[MEAN_MOLECULAR_WEIGTH];
  g_isoSoundSpeed = sqrt(cs02);
 #endif
#else
 #if EOS == ISOTHERMAL
  g_isoSoundSpeed = g_inputParam[ASPECT_RATIO] * sqrt(G_MU);
  cs = g_isoSoundSpeed;
 #endif
#endif

  // Calculate kinematic viscosity
 #if FIXED_ASPECT_RATIO == NO
  H = cs/OmegaK;
  g_nu = g_inputParam[ALPHA] * cs02 / sqrt(G_MU);
 #else
  H = g_inputParam[ASPECT_RATIO] * R;
  cs = H*OmegaK;
  g_nu = g_inputParam[ALPHA] * POW2(g_inputParam[ASPECT_RATIO]) * sqrt(G_MU);
 #endif

  us[RHO] = sigma; // /(sqrt(2.*CONST_PI)*H);
  Omega = OmegaK * sqrt((p_index+q_index)*POW2(H/R) + 1.0);

  // Define the gas velocities
  us[VX1] = us[VX2] = us[VX3] = 0.0;
  #if ROTATING_FRAME == YES
  g_OmegaZ = sqrt(G_MU);
  #endif
  us[iVPHI] = R*(Omega - g_OmegaZ);
  #ifdef FARGO
  us[FARGO_W] = R*(Omega - g_OmegaZ);
  us[iVPHI] = 0.0;
  #endif  

  // Define the pressure for the ideal EOS
  #if EOS == IDEAL
  g_gamma = 1.01;
  us[PRS] = us[RHO]*cs*cs;
  #endif
  us[TRC] = 0.0;
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{  
  Particle *part;
  particleNode *CurNode;

  int bin;
  bool size_cond, radius_cond;
  double r_part, phi_part, r_planet, phi_planet, dx, dy, dr2, ds2;

  double r_in = grid->xbeg_glob[IDIR]*g_inputParam[DAMPING_INNER];
  double r_out = grid->xend_glob[IDIR]*g_inputParam[DAMPING_OUTER];
  
  //cutoff radius (don't include particles in injection zone)
  double r_cutoff = r_out-(r_out-r_in)*g_inputParam[DUST_INJECTION_RATIO];

  // WARNING: the variables NBIN_DUST and N_BIN_DUST are very different
  // NBIN_DUST = number of bins ---> user-defined parameter in pluto.ini
  // N_BIN_DUST = 10 --------------> it's the index of the NBIN_DUST element in g_InputParam

  double local_dust_force[NBIN_DUST];
  double global_dust_force[NBIN_DUST];
  long int local_particles[NBIN_DUST];
  long int global_particles[NBIN_DUST];

  //initialize forces and counts
  
  for (bin = 0; bin < NBIN_DUST; bin++) 
  {
    local_dust_force[bin]=0.0;
    global_dust_force[bin]=0.0;
    local_particles[bin]=0;
    global_particles[bin]=0;
  }
  
  //planet polar coordinates
  phi_planet = atan2(g_nb.y[1],g_nb.x[1]);
  r_planet = sqrt(g_nb.x[1]*g_nb.x[1]+g_nb.y[1]*g_nb.y[1]);

  //potential smoothing length
  ds2=POW2(DUST_SMOOTHING_FACTOR*r_planet*pow(g_nb.m[1]/(3*g_nb.m[0]),1/3));

  //loop on particles
  PARTICLES_LOOP(CurNode, d->PHead){
    part = &(CurNode->p);
    r_part=part->radius*UNIT_LENGTH;
    //loop on particle size
    for (bin = 0; bin < NBIN_DUST; bin++) 
    {
      //one bin at a time (with rounding errors)
      size_cond = ( fabs(r_part-dust_size_array[bin]) < 1e-8 * fmax(fabs(r_part),(dust_size_array[bin])) ); 
      
      //only outisde injection zone (no artificial density evolution)
      radius_cond = (part->coord[IDIR]<r_cutoff); 

      if (size_cond && radius_cond) 
      {
        // if (bin==1){
        //   printf("r_part: %6e,  dust_bin: %6e",r_part,dust_size_array[bin]);
        // }
        //particle polar coordinates
        r_part = part->coord[IDIR];
        phi_part = part->coord[JDIR];
        
        //particle cartesian coordinates on rotating planet frame (x radial)
        dx=r_part*cos(phi_part-phi_planet)-r_planet;
        dy=r_part*sin(phi_part-phi_planet);

        //y force calculation

        dr2=dx*dx+dy*dy+ds2;
        local_dust_force[bin]+=dy/pow(dr2,1.5);
        local_particles[bin]++;
      }
    }
  }

  //merge the contributions from all processors
  MPI_Allreduce(local_dust_force, global_dust_force, NBIN_DUST, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(local_particles, global_particles, NBIN_DUST, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  
  // for (bin = 0; bin < NBIN_DUST; bin++) {
  //   print("bin: %i N particles: %lu; dust force %e\n", bin, global_particles[bin], global_dust_force[bin]);
  // }

  //output to file
  if (prank==0){
  
    char output_file[512];
    sprintf(output_file, "%s/dust_force.dat", RuntimeGet()->output_dir);

    //header
    if (g_stepNumber == 0)
    {
      FILE *fp = fopen(output_file, "w");
      fprintf(fp,"# Number of bins: %d", NBIN_DUST);
      fprintf(fp,"\n#\n");
      fprintf(fp,"# Size [cm]:   ");
      for(bin = 0; bin < NBIN_DUST; bin++) {
        fprintf(fp,"%14.6e", dust_size_array[bin]);
        }
      for(bin = 0; bin < NBIN_DUST; bin++) {
        fprintf(fp,"%15.6e", dust_size_array[bin]);
        }
      fprintf(fp,"\n#\n");
      fprintf(fp,"%-17s","# Time");
      fprintf(fp,"%-*s", 14*NBIN_DUST, "Normalized dust force");
      fprintf(fp,"Number of particles counted");
      fprintf(fp,"\n#\n");
      fclose(fp);
    }

    //iteration output
    FILE *fp = fopen(output_file, "a");
    fprintf(fp,"%-15.6e",g_time);
    for(bin = 0; bin < NBIN_DUST; bin++) {
      fprintf(fp,"% 14.6e",global_dust_force[bin]);
    } 
    for(bin = 0; bin < NBIN_DUST; bin++) {
      fprintf(fp,"%15ld",global_particles[bin]);
    }
    fprintf(fp,"\n");
    fclose(fp);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* 
 *
 *********************************************************************** */
{
  int i, j, k, nv;
  double *x1, *x2, *x3, R, OmegaK,v[256];
  double p_index = -g_inputParam[P_INDEX];
  double q_index = -g_inputParam[Q_INDEX];
  double unit_surface_density = CONST_Msun/(UNIT_LENGTH*UNIT_LENGTH);

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  double rho0, damping, H, sigma;
  double rbeg = grid->xbeg_glob[IDIR];
  double rend = grid->xend_glob[IDIR];
  double inner_damping_radius = rbeg*g_inputParam[DAMPING_INNER];
  double outer_damping_radius = rend*g_inputParam[DAMPING_OUTER];

  if (side == 0) { /* -- check solution inside the domain -- */
    DOM_LOOP(k,j,i){

     #if GEOMETRY == POLAR
      R = x1[i];
     #elif GEOMETRY == SPHERICAL
      R = x1[i]*sin(x2[j]);
     #endif

      if (R < inner_damping_radius) {
        damping = POW2((inner_damping_radius-R)/(inner_damping_radius-rbeg));
        OmegaK = sqrt(G_MU/R)/R;
        
        d->Vc[VX1][k][j][i] *= 1.0-damping*g_dt*OmegaK;
        
        sigma = (g_inputParam[SIGMA_0] / unit_surface_density) * pow(R, p_index);
       #if FIXED_ASPECT_RATIO == NO
        double temperature = g_inputParam[TEMPERATURE_0] * pow(R, q_index);
        double cs = sqrt(temperature/KELVIN/g_inputParam[MEAN_MOLECULAR_WEIGTH]);
        H = cs/OmegaK;
       #else
        H = g_inputParam[ASPECT_RATIO] * R;
       #endif
        rho0 = sigma; // /(sqrt(2.*CONST_PI)*H);
        d->Vc[RHO][k][j][i] = rho0 + (d->Vc[RHO][k][j][i]-rho0) * exp(-damping*g_dt*OmegaK);
      }

      if (R > outer_damping_radius) {
        damping = POW2((R-outer_damping_radius)/(rend-outer_damping_radius));
	OmegaK = sqrt(G_MU/R)/R;

        d->Vc[VX1][k][j][i] *= 1.0-damping*g_dt*OmegaK;
        
        sigma = (g_inputParam[SIGMA_0] / unit_surface_density) * pow(R, p_index);
       #if FIXED_ASPECT_RATIO == NO
        double temperature = g_inputParam[TEMPERATURE_0] * pow(R, q_index);
        double cs = sqrt(temperature/KELVIN/g_inputParam[MEAN_MOLECULAR_WEIGTH]);
        H = cs/OmegaK;
       #else
        H = g_inputParam[ASPECT_RATIO] * R;
       #endif
        rho0 = sigma; // /(sqrt(2.*CONST_PI)*H);
        d->Vc[RHO][k][j][i] = rho0 + (d->Vc[RHO][k][j][i]-rho0) * exp(-damping*g_dt*OmegaK);
      }
    }
  }

  if (side == X1_BEG){
    X1_BEG_LOOP(k,j,i){
      NVAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][j][2*IBEG - i - 1];
      d->Vc[VX1][k][j][i] *= -1.0;
      #if GEOMETRY == POLAR
      R = x1[i];
      #elif GEOMETRY == SPHERICAL
      R = x1[i]*sin(x2[j]);
      #endif
      OmegaK = sqrt(G_MU/R)/R;
      d->Vc[iVPHI][k][j][i] = R*(OmegaK - g_OmegaZ);
      #ifdef FARGO
      d->Vc[iVPHI][k][j][i] = 0.0;
      #endif
      #if DUST_FLUID == YES
      d->Vc[VX2_D][k][j][i] = d->Vc[iVPHI][k][j][i];
      #endif
    }
  }

  if (side == X1_END){
    X1_END_LOOP(k,j,i){
      NVAR_LOOP(nv)  d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IEND];
      #if GEOMETRY == POLAR
      R = x1[i];
      #elif GEOMETRY == SPHERICAL
      R = x1[i]*sin(x2[j]);
      d->Vc[iVTH][k][j][i] = 0.0;
      #endif
      d->Vc[iVR][k][j][i]  = 0.0;
      OmegaK = sqrt(G_MU/R)/R;
      d->Vc[iVPHI][k][j][i] = R*(OmegaK - g_OmegaZ);
      #ifdef FARGO
      d->Vc[iVPHI][k][j][i] = 0.0;
      #endif
      #if DUST_FLUID == YES
      d->Vc[VX2_D][k][j][i] = d->Vc[iVPHI][k][j][i];
      #endif
    }
  }
}

#if (BODY_FORCE & VECTOR)
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
}
#endif

