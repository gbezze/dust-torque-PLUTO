# Simulation list

### Numerical parameters

| *Name* |  | *INI file name* | *unit* | *default value* |
|---|---|---|---|---|
| Planet mass        |$M_p$| value in `planet.ini` |$M_\odot$|    |
| Gas density at 1 AU | $\Sigma_0$ | `SIGMA_0` | $\text{g}/\text{cm}^2$ | 1000 |
| Dust grain density | $\rho_\bullet$ | `DUST_DENSITY` | $\text{g}/\text{cm}^3$ | 2
| Number of particles|$N_P$|`Nparticles`|  | $3\cdot 10^6$
| Dust min size | $s_{min}$ | `DUST_SIZE_MIN` |  $\text{cm}$  |  0.02
| Dust max size | $s_{max}$ | `DUST_SIZE _MAX` |  $\text{cm}$  |  20
| Number of size bins | $N_s$ | `N_BIN_DUST` |  |  30
| Simulation duration | $t_{sim}$ | `tstop` | $2\pi/\Omega_p$ | 2000
| Analysis period | $\Delta t_a$ | `analysis` | $2\pi/\Omega_p$ | 0.01

### Simulation options

- **Dust diffusion:** Simulates effects of gas turbulence in dust drag by adding a random velocity term $\delta \vec v _t $ to drag. Set up in `definitions.h`, can be swithced on or off with the definition `PARTICLES_DUST_DIFFUSION`.

- **Potential smoothing:** Avoids singularities ad takes into account 3D geometry by adding an $r_s$ term in quadrature to the distance in the planet's gravitational potential. The options are:
    - *Hill radius*  $r_s = 0.6 \ r_H = 0.6\ r_p(M_p/3 M_\star)^{1/3}$ : Not realistic but done in the original paper
    - *Gas scale height*  $r_s = 0.7 \ H = 0.7 \ h \ r_p$ : More realistic for very well coupled particles
    - *Dust scale height*  $r_s = 0.7 \ H_d = 0.7 H \sqrt{\alpha/(\alpha+\mathcal{S})}$ : takes into account vertical settling of dust. Potentially more accurate

### Output folders

The all parameters that are not listed have the default value.

| ID | $M_P \ [M_\odot]$ | $r_s$ | diffusion |
|----|---:|----|:--:|
|A1 | $1\cdot10^{-6}$ | $0.6\ r_H$| `NO` |
|A2 | $3\cdot10^{-6}$ | $0.6\ r_H$| `NO` |
|A3 | $1\cdot10^{-5}$ | $0.6\ r_H$| `NO` |
|A4 | $3\cdot10^{-5}$ | $0.6\ r_H$| `NO` |


