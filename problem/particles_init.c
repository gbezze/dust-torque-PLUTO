/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Set particles initial conditions for the fluid-CR relative drift test.
 
 Initialize particles for the Fluid-particles relative drift test,
 section 4.3 of [MVBM18].
 
 \authors A. Mignone (mignone@ph.unito.it)\n
          B. Vaidya (bvaidya@unito.it)\n
 
 \date    March 20, 2018
  \b References: \n
   - [MVBM18] "A PARTICLE MODULE FOR THE PLUTO CODE: I - AN IMPLEMENTATION OF
               THE MHD-PIC EQUATIONS", Mignone etal.ApJS (2018)  [ Sec. 4.3 ]
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

double dust_radial_exp(double size){

  //empyrical power law exponent for stationary dust distribution
  //size in cm

  double x=size/2.38;
  double exponent = 2.387*log10(pow(x,1.041)+pow(x,-0.117));
  //printf("size: %e, exponent %e", size, exponent);
  return exponent;
}


double generate_random_radius(double size, double ri, double ro, double rb){

  //generate random radius from stationary distribution with rejection method

  int found = 0;
  int N_try = 0;
  double r, X, q;
  double pdf=0.;

  while ((found==0)&&(N_try<100)){

    r = ri+(rand()/(1.+RAND_MAX))*(ro-ri);
    X = (rand()/(1.+RAND_MAX)); //dust density normalized at 1 for r=rb

    q = dust_radial_exp(size);

    //stationary dust radial density 
    if (r<rb){
      pdf=pow(r/rb,q);
    }
    else{
        //injection zone distribution (0 at r_end and 1 at r_damping)
        pdf = pow(r/rb,-1*q) * (log10(r/ro)/log10(rb/ro));
    }
    
    //injection zone distribution linear cutoff (0 at r_end and 1 at r_damping)
    // double K = 1.;
    // if (r>rd){
    //   K=(ro-r)/(ro-rd);
    // }

    N_try++;
    
    if (X<pdf){
      found=1;
    }
  }

  return r;
}

double size_dust(int bin){
  return pow(SIZE_EXP_BASE, bin)*g_inputParam[DUST_SIZE_MIN];
}

/* ********************************************************************* */
void Particles_Init(Data *d, Grid *grid)
/*!
 *  Sets initial conditions on particles.
 *
 *  \param [in]    d       Pointer to the PLUTO data structure.
 *  \param [in]    grid    Pointer to the PLUTO grid structure.
 *
 *********************************************************************** */
{
  int i,j,k, np, dir, nc, bin;
  double p;
  int np_glob = RuntimeGet()->Nparticles_glob;
  int np_cell = RuntimeGet()->Nparticles_cell;
  Particle part;
  double xbeg[3], xend[3], R, H, gas_density;
  double outer_radius, inner_radius, injection_boundary, OmegaK;
  double unit_surface_density = CONST_Msun/(UNIT_LENGTH*UNIT_LENGTH);
  int n_count = 0;
  double n = -g_inputParam[P_INDEX];
  double cs, temperature, gas_sigma;

  for (dir = 0; dir < 3; dir++) {
    xbeg[dir] = grid->xbeg_glob[dir];
    xend[dir] = grid->xend_glob[dir];
  }
  
  outer_radius = xend[IDIR]*g_inputParam[DAMPING_OUTER];
  inner_radius = xbeg[IDIR]*g_inputParam[DAMPING_INNER];
  injection_boundary = outer_radius-(outer_radius-inner_radius)*g_inputParam[DUST_INJECTION_RATIO];

  for (bin = 0; bin < NBIN_DUST; bin++) {
    particle_count_outer[bin] = 0;
    if (bin >= INIT_BIN_DUST) {
      dust_size_array[bin-INIT_BIN_DUST] = size_dust(bin-INIT_BIN_DUST);
    }
  }

/* --------------------------------------------------------------
   1. Global initialization
   -------------------------------------------------------------- */
  if (np_glob > 0){

    for (np = 0; np < np_glob; np++){
      for (bin = INIT_BIN_DUST; bin < NBIN_DUST; bin++) {
        // Assign particles in NBIN_DUST log spaced size bins between DUST_SIZE_MIN and DUST_SIZE_MAX (in cm)
        if (np%(NBIN_DUST-INIT_BIN_DUST) == (bin-INIT_BIN_DUST)) {
          //part.radius = pow(2, (bin-INIT_BIN_DUST))*g_inputParam[DUST_SIZE_MIN]/UNIT_LENGTH;
          part.radius = dust_size_array[bin-INIT_BIN_DUST]/UNIT_LENGTH;
          // Intialize particles positions to match gas density
          // part.coord[IDIR] = pow((pow(xend[IDIR], n+1) - pow(xbeg[IDIR], n+1))*(rand()/(1.+RAND_MAX)) + pow(xbeg[IDIR], n+1), 1/(n+1));
          part.coord[IDIR] = generate_random_radius(part.radius*UNIT_LENGTH, inner_radius, outer_radius, injection_boundary);
          //printf("radius: %e \n",part.coord[IDIR]);
          if (part.coord[IDIR] >= injection_boundary) {
            particle_count_outer[bin] += 1;
          }
        } 
      }
      
      R = part.coord[IDIR];

      for (i=1;i<DIMENSIONS;++i){
        part.coord[i] = xbeg[i] + (xend[i]-xbeg[i])*(rand()/(1.+RAND_MAX));
      }

      // Assign keplerian velocity
      part.speed[IDIR]    = 0.0;
      part.speed[JDIR]    = sqrt(G_MU/R);
      part.speed[KDIR]    = 0.0;

      // Assign particle stopping time
      p = -g_inputParam[P_INDEX];
      gas_sigma = (g_inputParam[SIGMA_0] / unit_surface_density) * pow(R, p);
      OmegaK = sqrt(G_MU/R)/R;
      part.density = g_inputParam[DUST_DENSITY]/UNIT_DENSITY;
      part.tau_s = 0.5*CONST_PI*part.radius*part.density/(gas_sigma*OmegaK);

     #if PARTICLES_DUST_SNOWLINE == YES      
      part.density = g_inputParam[DUST_DENSITY]/UNIT_DENSITY;
      part.composition = 'I';
     #else
      part.density = g_inputParam[DUST_DENSITY]/UNIT_DENSITY;
      part.composition = 'R';
     #endif

      Particles_Insert (&part, d, PARTICLES_CREATE, grid);
    }
    if (prank==0){
      printf("PARTICLES GENERATED:\n");
      printf("| outer dust radius:  %e \n",outer_radius);
      printf("| inner dust radius:  %e \n",inner_radius);
      printf("| injection boundary: %e \n",injection_boundary);
      printf("| \n");
      for (bin = INIT_BIN_DUST; bin < NBIN_DUST; bin++){
        printf("| %li particles in injection zone for bin %i \n",particle_count_outer[bin],bin);
      }
      printf("=============================================\n");
    }

  }
  Particles_SetID(d->PHead);
}

/* ********************************************************************* */
void Particles_Inject(Data *d, Grid *grid)
/*!
 *  Sets user-defined boundary conditions on particles.
 *
 *  \param [in]  data    Pointer to the PLUTO data structure.
 *  \param [in]  grid    Pointer to the PLUTO grid structure.
 *********************************************************************** */
{

  int bin, dir, np, ipart;
  double p;
  double R, OmegaK, H, gas_density;
  Particle *part, pnew;
  double unit_surface_density = CONST_Msun/(UNIT_LENGTH*UNIT_LENGTH);
  //Gas domain radial boundaries
  double rend = grid->xend_glob[IDIR];
  double rbeg = grid->xbeg_glob[IDIR];

  //Dust domain radial boundaries
  double outer_radius = rend*g_inputParam[DAMPING_OUTER]; 
  double inner_radius = rbeg*g_inputParam[DAMPING_INNER];

  //Dust injection area boundary
  double injection_boundary = outer_radius - (outer_radius - inner_radius)*g_inputParam[DUST_INJECTION_RATIO];
  
  
  double gas_sigma;

  p = -g_inputParam[P_INDEX];

  Particles_Check_OuterFlux(d, grid);

  //min number of particles in outer bin
  // int min_particles_outer = particle_count_outer[1];

  // for (bin=INIT_BIN_DUST; bin<NBIN_DUST; bin++) {
  //   if (particle_count_outer[bin]>min_particles_outer) {
  //     min_particles_outer=particle_count_outer[bin];
  //   }
  // }

  for (bin=INIT_BIN_DUST; bin<NBIN_DUST; bin++) {

    if (global_count_outer[bin] < particle_count_outer[bin]) {

      for (ipart=0; ipart<(particle_count_outer[bin] - global_count_outer[bin]); ipart++) {

        pnew.radius = size_dust(bin)/UNIT_LENGTH;
        pnew.coord[IDIR] = injection_boundary + (outer_radius - injection_boundary)*(rand()/(1.+RAND_MAX));
        //pnew.coord[IDIR] = pow((pow(rend, p+1) - pow(outer_damping_radius, p+1))*(rand()/(1.+RAND_MAX)) + pow(outer_damping_radius, p+1), 1/(p+1));
        pnew.coord[JDIR] = 2.0*CONST_PI*(rand()/(1.+RAND_MAX));
        pnew.coord[KDIR] = 0.0;
        pnew.speed[IDIR] = 0.0;
        pnew.speed[JDIR] = sqrt(G_MU/pnew.coord[IDIR]);
        pnew.speed[KDIR] = 0.0;
        pnew.density = g_inputParam[DUST_DENSITY]/UNIT_DENSITY;;
        R = pnew.coord[IDIR];
        gas_sigma = (g_inputParam[SIGMA_0] / unit_surface_density) * pow(R, p);
        OmegaK = sqrt(G_MU/R)/R;
        pnew.tau_s = 0.5*CONST_PI*pnew.radius*pnew.density/(gas_sigma*OmegaK);
        pnew.composition = 'R';
        pnew.id = RuntimeGet()->Nparticles_glob + 1;

        Particles_Insert(&pnew, d, PARTICLES_CREATE, grid);
      }
    }
  }
  Particles_SetID(d->PHead);
}

/* ********************************************************************* */
void Particles_UserDefBoundary(Data *d, int side, Grid *grid)
/*
 *
 *********************************************************************** */
{
  int    i, dir;
  double xbeg[3], xend[3], r_inner, r_outer;
  particleNode *curr = d->PHead, *next;
  Particle *p;

  for (dir = 0; dir < 3; dir++) {
    xbeg[dir] = grid->xbeg_glob[dir];
    xend[dir] = grid->xend_glob[dir];
  }

  r_inner=xbeg[IDIR]*g_inputParam[DAMPING_INNER];
  r_outer=xend[IDIR]; //allow dust to go above dust injection area (by turbolent diffusion or kicks from planet)

  if (side == X1_BEG){
    dir = IDIR;
    while (curr != NULL){
      p    = &(curr->p);
      next = curr->next;
      if (p->coord[dir] < r_inner) {
        Particles_Destroy (curr,d);
      }
      curr = next;
    }
  }

  if (side == X1_END){
    dir = IDIR;
    while (curr != NULL){
      p    = &(curr->p);
      next = curr->next;
      if (p->coord[dir] >= r_outer) {
        Particles_Destroy (curr,d);
      }
      curr = next;
    }
  }

  if (side == X2_BEG){
    dir = JDIR;
    while (curr != NULL){
      p    = &(curr->p);
      next = curr->next;
      if (p->coord[dir] < xbeg[dir]) {
        p->coord[dir] += 2.*CONST_PI;
      }
      curr = next;
    }
  }

  if (side == X2_END){
    dir = JDIR;
    while (curr != NULL){
      p    = &(curr->p);
      next = curr->next;
      if (p->coord[dir] > xend[dir]) {
        p->coord[dir] -= 2.*CONST_PI;
      }
      curr = next;
    }
  }

}
