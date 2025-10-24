/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Tools required to define the Particle MPI struct and 
        Interpolating quantities from grid to particles.
 
 \authors   A. Mignone (mignone@ph.unito.it)\n
            B. Vaidya (bvaidya@unito.it)\n
  
 \date     Mar 17, 2021
 */
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

// Function to perform binary search
int binarySearch(double array[], int size, double target) {
  //search target in array, if found return the index

    int left = 0;
    int right = size - 1;

    double tolerance = target*1e-6;

    while (left <= right) {
        int middle = 0.5 * (left + right);

        // Check if the target is present at mid
        if ( (array[middle] < target + tolerance) && (array[middle] > target - tolerance) ) {
            return middle; // Target found
        }

        // If target is greater, ignore the left half
        if (array[middle] < target) {
            left = middle + 1;
        }
        // If target is smaller, ignore the right half
        else {
            right = middle - 1;
        }
    }
    return -1; // Target not found
}

/* ********************************************************************* */
void Particles_Check_OuterFlux(Data *d, Grid *grid)
/*!
 *
 *
 *********************************************************************** */
{
  int bin, il, im, ih;
  double r_out = grid->xend_glob[IDIR]*g_inputParam[DAMPING_OUTER];
  double r_in = grid->xbeg_glob[IDIR]*g_inputParam[DAMPING_INNER];
  double injection_boundary = r_out - (r_out-r_in) *g_inputParam[DUST_INJECTION_RATIO];
  Particle *part;
  particleNode *CurNode;

  for (bin = 0; bin < NBIN_DUST; bin++) {
    local_count_outer[bin] = 0;
  }

  PARTICLES_LOOP(CurNode, d->PHead){
    part = &(CurNode->p);
    if (part->coord[IDIR] > injection_boundary) {
      im = binarySearch(dust_size_array, NBIN_DUST-INIT_BIN_DUST, part->radius*UNIT_LENGTH);
      //index of the dust size in size array
      //print("BS: %le %le %d\n", part->radius*UNIT_LENGTH, dust_size_array[im], im);
      local_count_outer[im+INIT_BIN_DUST] += 1;
    }
  }

  MPI_Allreduce(local_count_outer, global_count_outer, NBIN_DUST, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

  // for (bin = INIT_BIN_DUST; bin < NBIN_DUST; bin++) {
  //  print("%d %d %d %d\n", bin, local_count_outer[bin], global_count_outer[bin], particle_count_outer[bin]);
  // }
}

/* ********************************************************************* */
long Particles_CheckAll (particleNode *PHead, int nbuf, Grid *grid)
/*!
 *  Count and return the number of particles inside a given region
 *  of the computational domain.
 *  By default (buf=0) the region is the active part of the domain (ghost
 *  zone excluded).
 *  This region may be enlarged by an additional buffer of nbuf zones in
 *  the ghost zones
 *
 *********************************************************************** */
{
  long int count=0, check;
  Particle     *p;
  particleNode *curr, *next;
  
  curr  = PHead;
  while (curr != NULL) {  /* Loop on particles */
    p    = &(curr->p);

    check = Particles_CheckSingle(p, nbuf, grid);
    if (check) count++;
    next = curr->next; 
    curr = next;
  }

  return count;
}

/* ********************************************************************* */
int Particles_CheckSingle(Particle *p, int nbuf, Grid *grid)
/*!
 * Check if the particle belongs to the local processor domain
 *
 *  \param [in]   p      pointer to Particle structure.
 *  \param [in]   nbuf   extra layer of ghost zones
 *                       nbuf = 0   check if particle is in active domain;
 *                       nbuf = 1   check if particle is in active zones + 1
 *                                  layer of ghost zones, and so forth;
 *                       
 *
 * \return TRUE if the particle is inside the specified region.
 *         FALSE otherwise.
 *********************************************************************** */
{
  int    dir, bbeg, bend, abeg, aend;
  int    cond  = 1;
  double xbeg, xend;
  
  DIM_LOOP(dir){
    abeg = grid->lbeg[dir];  /* Active domain starting index */
    aend = grid->lend[dir];  /* Active domain ending index */
    
    bbeg = abeg - nbuf;
    bend = aend + nbuf;

    xbeg = grid->xl[dir][bbeg];
    xend = grid->xr[dir][bend];

    cond *= (p->coord[dir] <= xend);
    cond *= (p->coord[dir] >  xbeg);   
  }

  return (cond != 0 ? 1:0);
}

/* ********************************************************************* */
double Particles_Interpolate(double ***V, double ***w, int *indx)
/*! 
 * Interpolate a grid quantity V to particle position x.
 *
 * \param [in]   V    a 3D array defined on the fluid grid->
 * \param [in]   w    a 3x3x3 array containing the weights
 * \param [in]  indx  a 3 element array giving the starting indices  (i,j,k). 
 *
 * \return The interpolated value of V at the desired location.
 *********************************************************************** */
{
  int    i, j, k;
  int    i1, j1, k1;
  double v = 0.0;
  
  i = indx[IDIR];
  j = indx[JDIR];
  k = indx[KDIR];
 
/* ----------------------------------------------------
    Interpolate 
   ---------------------------------------------------- */

  for (k1 = -INCLUDE_KDIR; k1 <= INCLUDE_KDIR; k1++){
  for (j1 = -INCLUDE_JDIR; j1 <= INCLUDE_JDIR; j1++){
  for (i1 = -INCLUDE_IDIR; i1 <= INCLUDE_IDIR; i1++){
    v += w[k1][j1][i1]*V[k+k1][j+j1][i+i1];
  }}}

  return v; 
}

/* ********************************************************************* */
void Particles_InterpolateArr (double ****V, int nfield, double ***w,
                               int *indx, double *v)
/*! 
 * Similar to ::Particles_Interpolate, but performs the interpolation on
 * n_field directly
 *
 * \param [in]  V        a 4D array containing 3D arrays defined on the grid
 * \param [in]  nfield   the number of arrays to be interpolated
 * \param [in]  w        a 3x3x3 array containing the weights
 * \param [in]  indx     a 3 element array giving the starting indices  (i,j,k). 
 * \param [out] v        a nfield element array containing interpolated values
 *
 *********************************************************************** */
{
  int    i, j, k, field;
  int    i1, j1, k1;
  
  i = indx[IDIR];
  j = indx[JDIR];
  k = indx[KDIR];
 
/* ----------------------------------------------------
    Interpolate 
   ---------------------------------------------------- */
  
  for (field = 0; field < nfield; field++){
    v[field] = 0.;
    for (k1 = -INCLUDE_KDIR; k1 <= INCLUDE_KDIR; k1++){
    for (j1 = -INCLUDE_JDIR; j1 <= INCLUDE_JDIR; j1++){
    for (i1 = -INCLUDE_IDIR; i1 <= INCLUDE_IDIR; i1++){
      v[field] += w[k1][j1][i1]*V[field][k+k1][j+j1][i+i1];
    }}}
  }
  
}

/* ********************************************************************* */
int Particles_LocateCell(double *xp, int *indx, Grid *grid) 
/*! 
 * Determine the index of the computational zone hosting the particle,
 * \f$ x_{i-\HALF} <= x_p < x_{i+\HALF} \f$.
 * The search extends to the entire computational domain, including
 * active and ghost zones.
 *    
 * This function works for both uniform and stretched grids and employs
 * a binary search algorithm.
 * 
 * \param [in]  xp    a 3-element array specifying the particle position
 * \param [out] indx  a 3-element array giving the cell coordinate (i,j,k)
 * \param [in]  grid  a pointer to an array of grid structures.
 *
 * Return 0 on success, 1 on error.
 *********************************************************************** */
{

  int i, dir, ngh;
  int l_ind, r_ind, m_ind;  
  double xL;

  indx[IDIR] = indx[JDIR] = indx[KDIR] = 0;

  // Binary search (any grid)

  DIM_LOOP(dir) {
    l_ind = 0;
    r_ind = grid->lend[dir] + grid->nghost[dir];
    while (l_ind != r_ind) {
      m_ind = l_ind + (r_ind - l_ind) / 2;
   
      if (xp[dir] < grid->xr[dir][m_ind]) {
        r_ind = m_ind;
      } else {
        l_ind = m_ind + 1;
      }
    }   
    indx[dir] = r_ind;

    if (indx[dir] < 0 || indx[dir] >= grid->np_tot[dir]){
      printLog ("! Particles_LocateCell(): particle outside the ");
      printLog ("computational domain\n");
      return 1;
    }  
  }
  return 0;

}


/* ********************************************************************* */
Particle *Particles_Select(particleNode *PHead, int id)
/*!
 *  Loop over particle and return the one with specified id.
 *
 * \param [in] PHead       pointer to the head node of the particle
 *                         linked list.
 * \param [in] id          the particle id
 *
 * Return the particle (if any), or NULL
 *********************************************************************** */
{
  particleNode *curNode;
  Particle     *p;

  PARTICLES_LOOP(curNode, PHead){
    p = &(curNode->p);
    if (p->id == id)  return p;
  }
  return NULL;
}

