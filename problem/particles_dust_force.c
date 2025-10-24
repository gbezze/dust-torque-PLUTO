/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute dust-related quantities (like force).

  \authors G. Picogna (picogna@usm.lmu.de)\n

  \b References
    - "MAGNETOHYDRODYNAMIC-PARTICLE-IN-CELL METHOD FOR COUPLING COSMIC RAYS 
       WITH A THERMAL PLASMA: APPLICATION TO NON-RELATIVISTIC SHOCKS"\n
       Bai et al., ApJ (2015) 809, 55

  \date   Nov 28, 2023
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

double calculate_disc_scale_height(double sound_speed, double radius)
{
  double H, OmegaK;

  OmegaK = sqrt(G_MU/radius)/radius;
  H = sound_speed/OmegaK;
  return H;
}

double calculate_density_midplane(double density, double disc_scale_height, double radius)
{
  double density_midplane;
 #if DIMENSIONS == 2
  density_midplane = density/(sqrt(2.*CONST_PI)*disc_scale_height);
 #elif DIMENSIONS == 3
  gas_density = density;
 #else
  #error
 #endif
  return density_midplane;
}

double calculate_gas_mean_free_path(double gas_density)
/* Calculate the mean free path for a 5:1 H2-He mixute.
   Lyra et al. 2009 eq. 19                              */
{
  double mu, sigma_mol, l;
  mu = 3.9e-24; // mean molecular weigth of a 5:1 H2-He mixture [g]
  sigma_mol = 2e-15; // cross section of molecular hydrogen [cm^2]
  l = (mu/sigma_mol/(gas_density*UNIT_DENSITY))/UNIT_LENGTH; // gas mean free path [cm]
  return l;
}

double calculate_drag_coefficient(double gas_density, double particle_radius, double mean_free_path, double Ma, double Kn, double c_s, double v_rel_abs)
/* 
   Calculate the drag coefficient interpolating between the Epstein and Stokes regimes.
   See Lyra et al. 2009 (https://www.aanda.org/articles/aa/pdf/2009/03/aa10797-08.pdf) eq. 12
*/
{
  double nu, Re, CdE, CdS, Cd;

  nu = sqrt(8./CONST_PI)*gas_density*c_s*mean_free_path/3.; // gas kinematic viscosity [g/cm/s]

  Re = 2.0*particle_radius*gas_density*v_rel_abs/nu;

  CdE = 2.0*sqrt(1.+128./9./CONST_PI/POW2(Ma));
  
  if (Re <= 500.0) {
    CdS = 24.0/Re + 3.6*pow(Re,-0.313);
  } else if (Re <= 1500.0) {
    CdS = 9.5e-5*pow(Re,1.397);
  } else {
    CdS = 2.61;
  }

  Cd = (9.0*Kn*Kn*CdE + CdS)/POW2(3.0*Kn+1.0);
  return Cd;
}

double calculate_sound_speed(Particle *p)
{
  double c_s;

 #if FIXED_ASPECT_RATIO == NO
  #if EOS == IDEAL
   c_s = sqrt(g_gamma*pressure/density);
  #elif EOS == ISOTHERMAL
   double q = -g_inputParam[Q_INDEX];
   double R = p->coord[0];
   c_s = g_isoSoundSpeed*pow(R, 0.5*q);
  #endif
 #else
  c_s = g_inputParam[ASPECT_RATIO] * sqrt(G_MU);
 #endif

  return c_s;
}

double calculate_temperature(double sound_speed)
{
  double T;
  T = POW2(sound_speed)*g_inputParam[MEAN_MOLECULAR_WEIGTH]*KELVIN;

  return T;
}

void Particles_Dust_UpdateSize(Particle *p, double gas_density)
{
  /*
    * Update the size of the dust particles based on the equilibrium vapor pressure
    * of water and the water vapor pressure.
    * See Garate et al. 2020 eqs. 33, 34
    * TODO: - extend to other species
  */
  double R, T, P_eq_water_0, A_water, P_eq_water;
  double mu_vap, h_g, P_water_vap, c_s;
  
  c_s = calculate_sound_speed(p);
  T = calculate_temperature(c_s);
  h_g = calculate_disc_scale_height(c_s, p->coord[0]);

  P_eq_water_0 = 1.14e13; // equilibrium vapor pressure of water at 0 K [dyn/cm^2] 
  A_water = 6062.0; // Clausius-Clapeyron constant [K]
  mu_vap = 18.0; // molecular weight of the water vapor [g/mol]
  
  // double A_CO, B_CO, P_eq_CO_0, P_eq_CO;
  // A_CO = 1030.0;
  // B_CO = 27.37;
  // P_eq_CO_0 = 1.0; // equilibrium vapor pressure of CO at 0 K [dyn/cm^2]
  // P_eq_CO = P_eq_CO_0 * exp(-A_CO/T + B_CO); // equilibrium vapor pressure of CO2 at 0 K [dyn/cm^2]

  P_eq_water = P_eq_water_0 * exp(-A_water/T);

  P_water_vap = gas_density * T / (sqrt(2.*CONST_PI) * h_g * mu_vap * KELVIN); // water vapor pressure [dyn/cm^2]

  if (P_water_vap < P_eq_water) { // when water vapour pressure is below this threshold the ice evaporates into vapour
    if (p->composition == 'I') {
      p->composition = 'R';
      p->radius *= 0.5;
      p->density *= 3.0;
    }
    //d->Vc[SIGMA_VAP-NFLX][k][j][i] += MIN(sqrt(2.*CONST_PI)*h_g*mu_vap*KELVIN*(P_eq_water-P_water_vap)/T, d->Vc[SIGMA_ICE-NFLX][k][j][i]);
  } else {
    if (p->composition == 'R') {
      p->composition = 'I';
      p->radius *= 2.0;
      p->density *= 0.3333333333333333;
    }
    //d->Vc[SIGMA_ICE-NFLX][k][j][i] += MIN(sqrt(2.*CONST_PI)*h_g*mu_vap*KELVIN*(P_water_vap-P_eq_water)/T, d->Vc[SIGMA_VAP-NFLX][k][j][i]);
  }
}

double Particles_Dust_StoppingTime(double v_rel_abs, double density, Particle *p)
{
  /* 
   Calculate the dust stopping (or friction) time based on the Epstein and Stokes regimes.
   See Lyra et al. 2009 (https://www.aanda.org/articles/aa/pdf/2009/03/aa10797-08.pdf) eq. 15
  */
  double c_s, mu, sigma_mol, l, f;
  double nu, Kn, Ma, Re;
  double CdE, CdS, Cd;
  double gas_density, h_g;

  c_s = calculate_sound_speed(p);
  h_g = calculate_disc_scale_height(c_s, p->coord[0]);

  gas_density = calculate_density_midplane(density, h_g, p->coord[0]);

  l = calculate_gas_mean_free_path(gas_density);

  Kn = 0.5*l/p->radius;
  Ma = v_rel_abs/c_s;

  Cd = calculate_drag_coefficient(gas_density, p->radius, l, Ma, Kn, c_s, v_rel_abs);

  return 4.0*l*p->density/(3.0*gas_density*Cd*c_s*Kn*Ma);
}

/* ******************************************************************************* */
void Particles_Dust_ComputeDrag(Particle *p, double *velocity, double density, double *drag_force)
/*
  * Calculate the drag acceleration onto the dust particles(g) and returns the
  * stopping time tstop.
  * TODO: - pure Stokes regime for reference
  *       - extend to other geometries
  *       - read values of radius a and mass m0 of the molecules from input file
  *       - calculate the geometric cross section sigma from a instead of using
  *         the fixed value
  ********************************************************************************* */
{
  double vrel[DIMENSIONS], vrelabs;

  #if ROTATING_FRAME == YES
  velocity[JDIR] += g_OmegaZ*p->coord[0];
  #endif
  DIM_EXPAND(vrel[IDIR] = velocity[IDIR] - p->speed[IDIR]; ,
              vrel[JDIR] = velocity[JDIR] - p->speed[JDIR]; ,
              vrel[KDIR] = velocity[KDIR] - p->speed[KDIR];
            )
  
  #if PARTICLES_DUST_STOPPING_TIME == VARIABLE
  vrelabs = DIM_EXPAND(vrel[0]*vrel[0], + vrel[1]*vrel[1], + vrel[2]*vrel[2]);
  vrelabs = sqrt(vrelabs);

  p->tau_s = Particles_Dust_StoppingTime(vrelabs, density, p);
  #endif

  DIM_EXPAND(drag_force[IDIR] = vrel[IDIR]/p->tau_s;,
              drag_force[JDIR] = vrel[JDIR]/p->tau_s;,
              drag_force[KDIR] = vrel[KDIR]/p->tau_s;
  )
}

/* *********************************************************************************** */
void transformCartesianPolar(double *cartesian, double phi, double *polar)
/*
 * Transform a vector from cartesian to polar coordinates.
 * ******************************************************************************* */
{
  DIM_EXPAND(polar[0]  = cartesian[0]*cos(phi) + cartesian[1]*sin(phi); ,
             polar[1]  =-cartesian[0]*sin(phi) + cartesian[1]*cos(phi); ,
             polar[2]  = cartesian[2];
  )
}

/* ******************************************************************************* */
/* ******************************************************************************* */

void Particles_Dust_ComputeGravity(Particle *p, double *grav_force)
/*
 * Calculate the gravitational acceleration onto the dust particles.
 * TODO: - extend to other geometries
 * *************************************************************************** */
{
  double r, r3, theta, phi, ds2;
  double grav_cart[DIMENSIONS], grav_polar[DIMENSIONS];
  int l;

  r = p->coord[0];
  phi = p->coord[1];
  r3 = POW3(r);

  DIM_EXPAND(double x = r * cos(phi);,
             double y = r * sin(phi);,
             double z = 0.0;)

  DIM_EXPAND(grav_cart[0] = -G_MU * x / r3;,
             grav_cart[1] = -G_MU * y / r3;,
             grav_cart[2] = -G_MU * z / r3;)

  for (l = CENTRAL_OBJECT; l < NB_N; l++)
  {
    double d, d3;
    DIM_EXPAND(
      double dx = x - g_nb.x[l];,
      double dy = y - g_nb.y[l];,
      double dz = z - g_nb.z[l];
    )

    d3 = pow(DIM_EXPAND(dx*dx, +dy*dy, +dz*dz)+ds2,1.5); // smoothing happens here
    
    DIM_EXPAND(
      grav_cart[0] += -g_nb.m[l]/g_nb.m[0] * G_MU * dx / d3;,
      grav_cart[1] += -g_nb.m[l]/g_nb.m[0] * G_MU * dy / d3;,
      grav_cart[2] += -g_nb.m[l]/g_nb.m[0] * G_MU * dz / d3;)
  }

  transformCartesianPolar(grav_cart, phi, grav_polar);

  for (l=0; l< DIMENSIONS; l++)
  {
      grav_force[l] = grav_polar[l];
  }
}

void Particles_Dust_Diffusion(Particle *p, double dt, Data *data, Grid *grid, double rho_gas)
/****************************************************************
 *
 * Add kicks in particle positions due to turbulence motions
 * (see Charnoz et al. 2011)
 *
 **************************************************************** */
{
    double pressure, Schmidt_number, Dd, c_s, H, R, tau, OmegaK;
    double vari, mean;
    double drho[DIMENSIONS], dD[DIMENSIONS];
    int i, j, k, l;

   #if GEOMETRY == SPHERICAL
    R = p->coord[0]*sin(p->coord[1]);
   #elif GEOMETRY == POLAR
    R = p->coord[0];
   #else
    print1("Geometry not implemented for the particle integrator/n")
    QUIT_PLUTO(1);
   #endif

    OmegaK = sqrt(G_MU/R)/R;

  #if FIXED_ASPECT_RATIO == NO
   #if EOS == IDEAL
    ParticlesInterpol(&pressure, pl, uu, grid, PRS);
    c_s = sqrt(g_gamma*pressure/rho_gas);
   #elif EOS == ISOTHERMAL
    c_s = g_isoSoundSpeed*pow(R, -0.5*g_inputParam[Q_INDEX]);
   #else
    print1("EOS not implemented for the particle integrator/n")
    QUIT_PLUTO(1);
   #endif
    H = c_s/OmegaK;
  #else
    H = g_inputParam[ASPECT_RATIO] * R;
    c_s = H*OmegaK;
  #endif

    tau = p->tau_s*OmegaK;
    Schmidt_number = POW2(1.0 + POW2(tau))/(1.0+4.0*tau); // Schmidt number
    Dd = g_inputParam[ALPHA]*c_s*H/Schmidt_number; // Dust diffusion coefficient

    i = p->cell[IDIR];
    j = p->cell[JDIR];
    k = p->cell[KDIR];

    DIM_EXPAND(
      drho[IDIR] = (data->Vc[RHO][k][j][i+1]-data->Vc[RHO][k][j][i])/grid->dx[IDIR][i];,
      drho[JDIR] = (data->Vc[RHO][k][j+1][i]-data->Vc[RHO][k][j][i])/grid->dx[JDIR][j];,
      drho[KDIR] = (data->Vc[RHO][k+1][j][i]-data->Vc[RHO][k][j][i])/grid->dx[KDIR][k];
    )

    // Compute the gradient of the diffusion coefficient
    // TODO: Ideal case
   #if EOS == ISOTHERMAL
    // In this case we know the analytical solution for the gradient of the Diffusion coefficient
    DIM_EXPAND(
      dD[IDIR] = (-g_inputParam[Q_INDEX]+1.5)*g_inputParam[ALPHA]/Schmidt_number*c_s*H*pow(R, -g_inputParam[Q_INDEX]+0.5);,
      dD[JDIR] = 0.;,
      dD[KDIR] = 0.
    )
   #else
    print1("EOS not implemented for the particle diffusion/n")
    QUIT_PLUTO(1);
   #endif

    for (l=0; l<DIMENSIONS; l++) {
      #if PARTICLES_VELOCITY_KICKS == YES
      vari = sqrt(2.0*Dd/dt + POW2(dD[l]));
      mean = Dd/rho_gas*drho[l] + dD[l];
      p->speed[l] += GaussianRandomNumber(mean, vari);
      #else // particle position kicks (default option)
      vari = sqrt(2.0*Dd*dt + POW2(dD[l]*dt));
      mean = (Dd/rho_gas*drho[l] + dD[l])*dt;
      p->coord[l] += GaussianRandomNumber(mean, vari);
      #endif
    }
}


/* ********************************************************************* */
void Particles_Dust_ComputeForce(Particle *p, double *v_gas, double rho_gas, double *force)
/*! 
 *********************************************************************** */
{
  double drag_force[DIMENSIONS], gravitational_force[DIMENSIONS];

  DIM_EXPAND(drag_force[IDIR] = 0.0;
             gravitational_force[IDIR] = 0.0;
             force[IDIR] = 0.0;,
             drag_force[JDIR] = 0.0;
             gravitational_force[JDIR] = 0.0;
             force[JDIR] = 0.0;,
             drag_force[KDIR] = 0.0;
             gravitational_force[KDIR] = 0.0;
             force[KDIR] = 0.0;)

  #if PARTICLES_DUST_DRAG == YES
  Particles_Dust_ComputeDrag(p, v_gas, rho_gas, drag_force);
  #endif
  #if PARTICLES_DUST_GRAVITY == YES
  Particles_Dust_ComputeGravity(p, gravitational_force);
  #endif

  DIM_EXPAND(force[IDIR] = drag_force[IDIR] + gravitational_force[IDIR];,
             force[JDIR] = drag_force[JDIR] + gravitational_force[JDIR];,
             force[KDIR] = drag_force[KDIR] + gravitational_force[KDIR];)
}
