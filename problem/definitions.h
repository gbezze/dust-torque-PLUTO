#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  GEOMETRY                       POLAR
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  USE_CMA                        YES
#define  USER_DEF_PARAMETERS            15

/* -- dust particles declarations -- */
#define  PARTICLES                      PARTICLES_DUST
#define  PARTICLES_DUST_TIME_STEPPING   SEMI_IMPLICIT
#define  PARTICLES_DUST_STOPPING_TIME   VARIABLE
#define  PARTICLES_DUST_DRAG            YES
#define  PARTICLES_DUST_GRAVITY         YES
#define  PARTICLES_DUST_DIFFUSION       YES
#define  PARTICLES_DUST_SNOWLINE        NO
#define  PARTICLES_VELOCITY_KICKS       NO
#define  PARTICLES_CHECK_WEIGHTS        YES
#define  PARTICLES_SHAPE                22

/* -- nbody declarations -- */
#define  NBODY_SYS                      YES
#define  CENTRAL_OBJECT                 STAR
#define  NO_OF_PLANETS                  1
#define  PLANET_FORMAT                  ORBIT
#define  CO_FEELS_DISK                  YES
#define  INDIRECT_TERM                  YES

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            ISOTHERMAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      EXPLICIT
#define  ROTATING_FRAME                 NO
#define  FIXED_ASPECT_RATIO             YES

/* -- user-defined parameters (labels) -- */

#define  M_CO                           0
#define  SIGMA_0                        1
#define  P_INDEX                        2
#define  Q_INDEX                        3
#define  ALPHA                          4
#define  ASPECT_RATIO                   5
#define  MEAN_MOLECULAR_WEIGTH          6
#define  DAMPING_INNER                  7
#define  DAMPING_OUTER                  8
#define  DUST_DENSITY                   9
#define  N_BIN_DUST                     10
#define  DUST_SIZE_MIN                  11
#define  DUST_SIZE_MAX                  12
#define  DUST_SMOOTHING_FACTOR          13
#define  DUST_INJECTION_RATIO           14

/* [Beg] user-defined constants (do not change this line) */

#define  INTERNAL_BOUNDARY              YES
#define  LIMITER                        VANLEER_LIM
#define  UNIT_LENGTH                    (CONST_au)
#define  UNIT_DENSITY                   (CONST_Msun/(UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH))
#define  UNIT_VELOCITY                  (sqrt(CONST_G*g_inputParam[M_CO]*CONST_Msun/UNIT_LENGTH)/(2.*CONST_PI))
#define  DEBUG                          NO
#define  MULTIPLE_LOG_FILES             NO
#define  CHOMBO_REF_VAR                 RHO
#define  SIZE_EXP_BASE                  pow(g_inputParam[DUST_SIZE_MAX]/g_inputParam[DUST_SIZE_MIN],1/(g_inputParam[N_BIN_DUST]-1))

/* [End] user-defined constants (do not change this line) */