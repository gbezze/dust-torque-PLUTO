# Simulation list

### Numerical parameters

| *Name* |  | *INI file name* | *unit* | *default value* |
|---|---|---|---|---|
| Planet mass        |$M_p$| value in `planet.ini` |$M_\odot$|    |
| Gas density at 1 AU | $\Sigma_0$ | `SIGMA_0` | $\text{g}/\text{cm}^2$ | 1000 |
| Dust grain density | $\rho_\bullet$ | `DUST_DENSITY` | $\text{g}/\text{cm}^3$ | 2
| Number of particles|$N_P$|`Nparticles`|  | $5 \cdot 10^5$
| Dust min size | $s_{min}$ | `DUST_SIZE_MIN` |  $\text{cm}$  |  0.2
| Dust max size | $s_{max}$ | `DUST_SIZE _MAX` |  $\text{cm}$  |  20
| Number of size bins | $N_s$ | `N_BIN_DUST` |  |  20
| Simulation duration | $t_{sim}$ | `tstop` | $2\pi/\Omega_p$ | 750
| Analysis period | $\Delta t$ | `analysis` | $2\pi/\Omega_p$ | 0.1
| Accretion timescale | $t_a$ | `DUST_ACCRETION_TIMESCALE` | $2\pi/\Omega_p$ | 1.0
### Simulation domain

The radial domains of the simulation are roughly proportional to the corotation half-width $\Delta r_c$:

|$M_p/M_\star$| $\Delta r_c/r_p$|$r_{ig}$ [AU] |$r_{og}$ [AU] |$\Delta_d$|
|--|--|--|--|--|
|$3 \cdot 10^{-5}$| $2.57 \cdot 10^{-2}$ | 0.75 | 1.25 | 0.05 |
|$1 \cdot 10^{-5}$| $1.48 \cdot 10^{-2}$ | 0.85 | 1.15 | 0.03 |
|$3 \cdot 10^{-6}$| $8.13 \cdot 10^{-3}$ | 0.92 | 1.08 | 0.015 |
|$1 \cdot 10^{-6}$| $4.7 \cdot 10^{-3}$ | 0.95 | 1.05 | 0.01 |

### Simulation options

- **Dust diffusion:** Simulates effects of gas turbulence in dust drag by adding a random velocity term $\delta \vec v _t $ to drag. Set up in `definitions.h`, can be swithced on or off with the definition `PARTICLES_DUST_DIFFUSION`.

- **Potential smoothing:** Avoids singularities ad takes into account 3D geometry by adding an $r_s$ term in quadrature to the distance in the planet's gravitational potential. The options are:
    - *Hill radius*  $r_s = \ r_H = \ r_p(M_p/3 M_\star)^{1/3}$ : Not realistic but done in the original paper
    - *Gas scale height*  $r_s = 0.7 \ H = 0.7 \ h \ r_p$ : More realistic for very well coupled particles
    - *Dust scale height*  $r_s = 0.7 \ H_d = 0.7 H \sqrt{\alpha/(\alpha+\mathcal{S})}$ : takes into account vertical settling of dust. Potentially more accurate

### Output folders

The all simulation parameters are encoded in the folder name:
`problem/out_[Smoothing][Planet mass]_[modifiers]/`
* `[Smoothing]` refers to the smoothing radius used:
    * `A`: Hill radius
    * `B`: Gas scale height
    * `C`: Dust scale height
    * `D`: No smoothing ($r_s=0.01r_H$)

* `[Planet mass]` is the planet mass in Earth Masses $M_\oplus$:
    $$M_\oplus/M_\odot = 3.33 \cdot 10^{-6}$$
    |$M_p/M_\star$|$M_\oplus$|
    |--:|--|
    |$3\cdot10^{-5}$|$9$|
    |$2\cdot10^{-5}$|$6$|
    |$1.5\cdot10^{-5}$|$4.5$|
    |$1\cdot10^{-5}$|$3$|
    |$7\cdot10^{-6}$|$2.1$|
    |$5\cdot10^{-6}$|$1.5$|
    |$3.33\cdot10^{-6}$|$1$|
    |$2\cdot 10^{-6}$|$0.6$|
    |$1.5 \cdot 10^{-6}$|$0.45$|
    |$1 \cdot {10^-6}$|$0.3$|

* `[modifiers]` is any parameter that differs form the standard ones:
    * `nd`: no diffusion
    * `na`: no accretion
    * `ta01`: short accretion time $t_a=0.1 \cdot2\pi \Omega_p^{-1}$
    * `ta2`: long accretion time $t_a=2 \cdot2\pi \Omega_p^{-1}$
    * `ld`: linear drag (Stokes regimes deactivated)