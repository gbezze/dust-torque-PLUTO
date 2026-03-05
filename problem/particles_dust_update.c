/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Update Dust particles.
 
  This file contains routines to update dust particles.
  
  \authors G. Picogna (picogna@usm.lmu.de)\n

  \b References
    - "PARTICLE CONCENTRATION AT PLANET-INDUCED GAP EDGES AND VORTICES. 
       I. INVISCID THREE-DIMENSIONAL HYDRO DISKS"\n
       Zhu et al., ApJ (2014) 785, 122

  \date   Nov 09, 2023
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

// the only change from the original file is that the planet position is saved to the g_inputparam array (line 157)

/* ********************************************************************* */
void Particles_Dust_Update(Data *data, timeStep *Dts, double dt, Grid *grid)
/*!
 * Advance particle by a step dt.
 * 
 * \param [in,out]  data  Data structure (contains particles list)
 * \param [in,out]  Dts   timeStep structure
 * \param [in]      dt    time increment
 * \param [in]      grid  pointer to Grid structure
 *********************************************************************** */
{
  int    dir;
  int    kcycle;
  double v_gas[DIMENSIONS], rho_gas, OmegaK;
  double dt_half, inv_dtnew, dts, dt_part;
  double centrifugal_force, force[DIMENSIONS], force_prev[DIMENSIONS];
  double coord_old[DIMENSIONS], speed_old[DIMENSIONS];
  double ts_prev = 0.;
  double inv_dt = 1.e-20;
  double centrifugal_force_old;
  static double ***w;
  static int first_call = 1;

  particleNode *CurNode;
  Particle *p;

  //DEBUG_FUNC_BEG ("Dust_Update");

  if (g_time < Dts->particles_tstart) return;

  #if SHOW_TIMING
  clock_t clock_beg = clock(), clock0;
  #endif

  #if GEOMETRY != POLAR
    #error Dust Particles work in POLAR geometry only
  #endif

  #if PARTICLES_DUST_DIFFUSION == YES
  if (first_call){
    RandomSeed(time(NULL),0);
    first_call = 0;
  }
  #endif

  Boundary (data, ALL_DIR, grid);

/* --------------------------------------------------------
   0. Allocate memory
   -------------------------------------------------------- */
  
  if (w == NULL){
    w      = ARRAY_BOX (-1, 1, -1, 1, -1, 1, double);
  }

/* --------------------------------------------------------
   1. Initialize sub-cycling
   -------------------------------------------------------- */

  dt      = dt/(double)Dts->Nsub_particles;
  dt_half = 0.5*dt;

/* --------------------------------------------------------
   2. Start sub-cycling
   -------------------------------------------------------- */

  Dts->invDt_particles = 1.e-38;

  for (kcycle = 0; kcycle < Dts->Nsub_particles; kcycle++){

    PARTICLES_LOOP(CurNode, data->PHead){

      p = &(CurNode->p);

    /* -- Save positions and velocities at the beginning of the step -- */
      for (dir = 0; dir < DIMENSIONS; dir++) {
        coord_old[dir] = p->coord[dir];
        speed_old[dir] = p->speed[dir];
        force_prev[dir] = 0.;
        force[dir] = 0.;
      }

    #if (PARTICLES_DUST_TIME_STEPPING == FULLY_IMPLICIT)
        /* -- Compute weights and indices at x^n ------------------------- */
        Particles_GetWeights(p, p->cell, w, grid);

        /* -- Interpolate fluid density and velocity at x^n -------------- */
        rho_gas = Particles_Interpolate(data->Vc[RHO], w, p->cell);
        DIM_EXPAND(
          v_gas[IDIR] = Particles_Interpolate(data->Vc[VX1], w, p->cell);,
          v_gas[JDIR] = Particles_Interpolate(data->Vc[VX2], w, p->cell);,
          v_gas[KDIR] = Particles_Interpolate(data->Vc[VX3], w, p->cell);)

        /* -- Compute forces acting on the dust particles at x^n --------- */
        Particles_Dust_ComputeForce(p, v_gas, rho_gas, force_prev);
        ts_prev = p->tau_s;

        /* -- Predict positions at x^{n+1} ------------------------------- */
        DIM_EXPAND(
          p->coord[IDIR] += p->speed[IDIR]*dt;,
          p->coord[JDIR] += p->speed[JDIR]/p->coord[IDIR]*dt;,
          p->coord[KDIR] += p->speed[KDIR]*dt;)
    #elif (PARTICLES_DUST_TIME_STEPPING == SEMI_IMPLICIT)
        /* -- Half drift of positions ------------------------------------ */
        DIM_EXPAND(
          p->coord[IDIR] += speed_old[0]*dt_half; ,
          p->coord[JDIR] += 0.5*(speed_old[1]/coord_old[0] + p->speed[1]/p->coord[0])*dt_half; ,
          p->coord[KDIR] += speed_old[2]*dt_half;)
    #else
        #error Particles_Dust_Update(): unknown integrator
    #endif

    #if PARTICLES_DUST_DIFFUSION == YES
      OmegaK = sqrt(G_MU/p->coord[IDIR])/p->coord[IDIR];
      // Assuming tau_eddy = 1/OmegaK and the probability of having a random kick dt/tau_eddy
      if (RandomNumber(0,1) < OmegaK * dt){
        /* -- Apply dust diffusion -- */
        Particles_GetWeights(p, p->cell, w, grid);
        rho_gas = Particles_Interpolate(data->Vc[RHO], w, p->cell);
        Particles_Dust_Diffusion(p, dt, data, grid, rho_gas);
      }
    #endif

      /* -- Compute weights and indices at new position ------------------ */
      Particles_GetWeights(p, p->cell, w, grid);

      /* -- Interpolate fluid density and velocity at new position ------- */
      rho_gas = Particles_Interpolate(data->Vc[RHO], w, p->cell);
      DIM_EXPAND(
        v_gas[IDIR] = Particles_Interpolate(data->Vc[VX1], w, p->cell);,
        v_gas[JDIR] = Particles_Interpolate(data->Vc[VX2], w, p->cell);,
        v_gas[KDIR] = Particles_Interpolate(data->Vc[VX3], w, p->cell);)

      /* -- Compute forces acting on the dust particles at new position -- */
      Particles_Dust_ComputeForce(p, v_gas, rho_gas, force);
      
      /* -- Save planet position to g_inputParam so that it can be used at the same timestep -- */
      g_inputParam[X_PLANET]=g_nb.x[1];
      g_inputParam[Y_PLANET]=g_nb.y[1];
      
    
      /* -- Calculate velocity update -- */
    #if (PARTICLES_DUST_TIME_STEPPING == FULLY_IMPLICIT)
        dts = 1.+dt_half*(p->tau_s + ts_prev + dt)/(ts_prev*p->tau_s);
       #if DIMENSIONS == 3
        p->speed[KDIR] += dt_half/dts*(force[KDIR] + force_prev[KDIR]*(1.+dt/ts_prev));
       #endif
        centrifugal_force_old = POW2(p->speed[JDIR])/coord_old[IDIR];
        p->speed[JDIR] += dt_half/dts*(force[JDIR] + force_prev[JDIR] * (1.+dt/ts_prev));
        centrifugal_force = POW2(p->speed[JDIR])/p->coord[IDIR];
        print("Before shift: %le %le %le %le %le %le\n", p->speed[IDIR]);
        p->speed[IDIR] += dt_half/dts*(force[IDIR] + centrifugal_force_old + (force_prev[IDIR] + centrifugal_force)*(1.+dt/ts_prev));
        print("After shift: %le %le %le %le %le %le %le\n", p->speed[IDIR], dt_half/dts, dt/ts_prev, force[IDIR], centrifugal_force_old, centrifugal_force, force_prev[IDIR]);
    #elif (PARTICLES_DUST_TIME_STEPPING == SEMI_IMPLICIT)
        dts = 1.+dt_half/p->tau_s;
       #if DIMENSIONS == 3
        p->speed[KDIR] = p->speed_old[KDIR] + force[KDIR]*dt/dts;
       #endif
        p->speed[JDIR] = speed_old[JDIR] + force[JDIR]*dt/dts;
        centrifugal_force = 0.5*(POW2(speed_old[JDIR])/coord_old[IDIR] + POW2(p->speed[JDIR])/p->coord[IDIR]);
        force[IDIR] += centrifugal_force;
        p->speed[IDIR] = speed_old[IDIR] + force[IDIR]*dt/dts;
    #else
        #error Particles_Dust_Update(): unknown integrator
    #endif

      /* -- Update spatial coordinate to obtain x^{n+1} -- */

    #if (PARTICLES_DUST_TIME_STEPPING == FULLY_IMPLICIT)
        p->coord[IDIR] = coord_old[IDIR] + dt_half*(speed_old[IDIR] + p->speed[IDIR]);
        p->coord[JDIR] = coord_old[JDIR] + dt_half*(speed_old[JDIR] + p->speed[JDIR]);
        // printLog ("after step: %le %le %le %le\n", p->coord[IDIR], dt_half, speed_old[IDIR], p->speed[IDIR]);
       #if DIMENSIONS == 3
        p->coord[KDIR] = coord_old[KDIR] + dt_half*(speed_old[KDIR] + p->speed[KDIR]);
       #endif
       // printLog ("After fully-implicit integrator... %le %le %le %le\n", coord_old[IDIR], dt_half, speed_old[IDIR], p->speed[IDIR]);
    #elif (PARTICLES_DUST_TIME_STEPPING == SEMI_IMPLICIT)
        p->coord[IDIR] += dt_half*p->speed[IDIR];
        p->coord[JDIR] += dt_half*p->speed[JDIR];
       #if DIMENSIONS == 3
        p->coord[KDIR] += dt_half*p->speed[KDIR];
       #endif
       // printLog ("After semi-implicit integrator... %le %le %le\n", p->tau_s, dt, p->coord[0]);
    #else
      #error Particles_Dust_Update(): unknown integrator
    #endif // SEMI_IMPLICIT

    /* ---------------------------------------------------------
        I1. Compute time step restriction 
       --------------------------------------------------------- */

      inv_dtnew = 0.2/p->tau_s;
      inv_dt = MAX(inv_dt, inv_dtnew);

    /* -- J. Check that particle has not travelled more than one cell --  */

      int checkp = 1;
      // checkp = Particles_CheckSingle(p, grid->nghost[IDIR], grid);
      if (checkp == 0){
        printLog ("! Particles_Dust_Update(): particle id %d outside domain\n",
                   p->id);
        printLog ("Active Domain: [%12.6e, %12.6e] [%12.6e, %12.6e]\n",
                   grid->xl[IDIR][IBEG], grid->xr[IDIR][IEND],
                   grid->xl[JDIR][JBEG], grid->xr[JDIR][JEND]);
        printLog ("Total  Domain: [%12.6e, %12.6e] [%12.6e, %12.6e]\n",
                   grid->xl[IDIR][0], grid->xr[IDIR][NX1_TOT-1],
                   grid->xl[JDIR][0], grid->xr[JDIR][NX2_TOT-1]);

        printLog ("Old coord: \n");
        ShowVector (coord_old, 3);
        printLog ("New coord: \n");
        ShowVector (p->coord, 3);
        double dxp[3];
        dxp[IDIR] = fabs(coord_old[IDIR] - p->coord[IDIR])/grid->dx[IDIR][IBEG];
        dxp[JDIR] = fabs(coord_old[IDIR] - p->coord[IDIR])/grid->dx[JDIR][JBEG];
        dxp[KDIR] = fabs(coord_old[IDIR] - p->coord[IDIR])/grid->dx[KDIR][KBEG];
        printLog ("|x(new) - x(old)|/dx: \n");
        ShowVector (dxp, 3);

        QUIT_PLUTO(1);
      }
      if (PARTICLES_DUST_SNOWLINE){
        Particles_Dust_UpdateSize(p, rho_gas);
      }
    }  /* End loop on particles */

  /* ----------------------------------------------------
     3e. Set boundary condition after deposition at
         x^{n+1/2} has been done
     ---------------------------------------------------- */

    Particles_Boundary(data, grid);
    Particles_BoundaryExchange(data, grid);

    #if SHOW_TIMING
    clock0 = clock();
    #endif
  }  /* End loop on sub-cycles */

  Dts->invDt_particles = inv_dt;

#if SHOW_TIMING
{
  double dclock_tot;

  Dts->clock_particles = (double)(clock() - clock_beg)/CLOCKS_PER_SEC;
  dclock_tot = (double)(clock() - clock_beg)/CLOCKS_PER_SEC;
  
  printLog ("  Total: %f, [%8.3e per particle]\n",dclock_tot,
                                               dclock_tot/p_nparticles);
}
#endif
  //DEBUG_FUNC_END ("Dust_Update");
}
