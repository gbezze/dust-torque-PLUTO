#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from particles_read import NStepStr
from scipy.stats import norm
from scipy.optimize import curve_fit
from pathlib import Path
from matplotlib import font_manager as fm, cycler
import os

global output_dir

output_dir = r'./problem/out_stat/'

def ReadPartData(ns):

    nstepstr = NStepStr(ns)
    fname="particles."+nstepstr+".partdbl"
    pathname=output_dir+fname
    h_lines = 0 
    val_dict = {}

	#READ HEADER. 
    with open(pathname, "rb") as f:
        for line in f:
            if h_lines < 13:
                if h_lines > 0 and h_lines < 13:
                    val_dict.update({line.split()[1].decode('utf8'):[i.decode('utf8') for i in line.split()[2:]]})
                h_lines += 1

    fdata = open(pathname, "rb")
 
    cnt = 0
    while (cnt < h_lines):
        fdata.readline()
        cnt += 1
        
    data = fdata.read()
    fdata.close()

    dt = np.dtype({'names':val_dict['field_names'], 'formats':['('+i+',)<d' for i in val_dict['field_dim']]})

    val = np.frombuffer(data, dtype=dt)

    for i in range(len(val_dict['field_names'])):
        name = val_dict['field_names'][i]
        val_dict.update({name:val[name]})

    return val_dict

def nice_plots(
    font_dir: str = r"/home/bezze/ebgaramond",
    palette: list[str] = None,
    font_size: int = 14,
    use_math: bool = True,
    ):
    """
    Configure Matplotlib to use EB Garamond with an optional custom color palette and math setup.

    Parameters
    ----------
    font_dir : str
        Directory containing EB Garamond .ttf files (Regular, Italic, Bold, etc.).
    palette : list of str
        Custom list of hex color codes for plotting.
    font_size : int
        Base font size for text.
    use_math : bool
        If True, sets mathtext to use Garamond italic / roman variants.
    """

    font_dir = Path(font_dir)
    if not font_dir.exists():
        raise FileNotFoundError(f"Font directory not found: {font_dir}")

    # --- Register fonts ---
    for font_file in font_dir.glob("EBGaramond*.ttf"):
        fm.fontManager.addfont(str(font_file))

    # --- Set base font ---
    plt.rcParams["font.family"] = "EB Garamond"
    plt.rcParams["font.size"] = font_size
    plt.rcParams["axes.titleweight"] = "bold"

    # --- Color palette ---
    if palette is None:
        palette = ["firebrick", "mediumseagreen", "darkorchid", "gold", "deepskyblue"]
    plt.rcParams["axes.prop_cycle"] = cycler(color=palette)

    # colormap
    colors = ["#AEEB52", "#DB8640", "#AD2B6A"]
    custom_cmap = mcolors.LinearSegmentedColormap.from_list("sanmartino", colors, N=256)
    try:
        plt.register_cmap(name="sanmartino", cmap=custom_cmap)
    except:
        pass

    plt.set_cmap(custom_cmap)
    # --- Math font settings ---
    if use_math:
        plt.rcParams.update({
            "mathtext.fontset": "custom",
            "mathtext.rm": "EB Garamond",
            "mathtext.it": "EB Garamond:italic",
            "mathtext.bf": "EB Garamond:bold",
            "mathtext.cal": "cmsy10",   # restore standard mathcal font
        })

    print("✅ EB Garamond style configured for Matplotlib.")

    return custom_cmap

def positions(step):

    # import dust data
    P = ReadPartData(step)

    r_part = np.array(P['x1']).flatten()
    phi_part = np.array(P['x2']).flatten()
    tau_part = np.array(P['tau_s']).flatten()
    size_part = np.array(P['radius']).flatten()

    time = float(P['time'][0])

    #import planet position

    time_index =step*2+10

    planet_data_path = os.path.join(output_dir, 'nbody_coordinates.dat')

    with open(planet_data_path,'r') as f:

        lines = f.readlines()
        x_planet = float(lines[time_index].split()[2])
        y_planet = float(lines[time_index].split()[3])
        t_planet = float(lines[time_index].split()[1])

    x_part = r_part * np.cos(phi_part) 
    y_part = r_part * np.sin(phi_part)

    truncate = lambda x, n: np.where(x==0, 0, np.trunc(x/10**(np.floor(np.log10(np.abs(x)))-n+1))*10**(np.floor(np.log10(np.abs(x)))-n+1))
    
    size_part = truncate(size_part,7)
    size_bins = np.sort(np.unique(size_part))

    # calculate individual particle tangential force on planet

    return x_part, y_part, size_part, size_bins , time, x_planet, y_planet, t_planet, tau_part

def read_force():
    
    force_data_path = os.path.join(output_dir, 'dust_force.dat')

    with open(force_data_path,'r') as f:

        lines = f.readlines()
        N_bins = int(lines[0].split(": ")[1].strip())
        N_times = len(lines)-6

        time=np.zeros(N_times)
        force=np.zeros((N_times,N_bins))
        counts=np.zeros((N_times,N_bins))

        bins=np.zeros(N_bins)
        bins = np.array(lines[2].split()[3:N_bins+3],dtype=float)

        for i in range(N_times):
            time[i]=lines[i+6].split()[0]
            force[i,:]=lines[i+6].split()[1:N_bins+1]
            counts[i,:]=lines[i+6].split()[N_bins+1:2*N_bins+1]

    return time, force, counts, N_bins, bins

def radial_distribution(S,step_start,step_end, plotting=True):

    from scipy.optimize import curve_fit

    N_steps = step_end - step_start
    radii=np.zeros(N_steps)

    for i in range(N_steps):

        step= step_start +i
        x, y, size_part, size_bins, time, x_planet, y_planet, t_planet, tau_part = positions (step)

        size_index = (size_part==size_bins[S])
        r= np.sqrt(x[size_index]**2 +y[size_index]**2)
        radii=np.append(radii,r)

    # fit power law, exclude damping zones

    # rmax = 1.5
    # rmin = 0.7
    # r_outer = 0.9 * rmax
    # r_inner = 1.1 * rmin

    rmax = 1.8
    rmin = 0.5
    r_outer = 0.9 * rmax
    r_inner = 1.1 * rmin
    r_outer= r_outer-(r_outer-r_inner)*0.1

    if plotting:
        plt.figure()
        counts, bins, _ = plt.hist(radii[(radii<r_outer) & (radii>r_inner)], bins=100, density=True)
    else:
        counts, bins = np.histogram(radii[(radii<r_outer) & (radii>r_inner)], bins=100, density=True)

    bin_centers = 0.5 * (bins[1:] + bins[:-1])

    def power_law(x, a, b):
        return a * x**b
    
    popt, pcov = curve_fit(power_law, bin_centers, counts, maxfev=10000)

    n0=popt[0]
    exponent=popt[1]
    exp_err=np.sqrt(pcov[1,1])

    if plotting:

        plt.plot(bin_centers, n0*bin_centers**exponent, 'r-', label=f'Fit: n0={n0:.2e}, exponent={exponent:.2f}')
        plt.legend()
        plt.show()
    
    return exponent, exp_err

#%% q_slope calculation

cmap=nice_plots()

time, force, counts, N_bins, bin_sizes = read_force()

exps = np.zeros(N_bins)
errs = np.zeros(N_bins)
for i in range(N_bins):
    exps[i], errs[i] = radial_distribution(i,17,20,plotting=False)

def dust_radial_exponent(x,x0,a,b,n):
    #for model fitting
    
    u=x/x0
    s1=u**(a*n)
    s2=u**(b*n)

    exponent= np.log10(s1+s2)/n

    return exponent

def dust_radial_exponent_p(x):
    #fitted model

    u=x/3.369
    return 2.246*np.log10(u**1.000+u**-0.124)

popt, pcov = curve_fit(dust_radial_exponent,bin_sizes,exps,[2,-0.2,1,0.4],errs)

print("")
print("OPTIMAL PARAMETERS for q(S)= K log10( (s/s0)^a + (s/s0)^b )")
print("")
print(" K  =  "+str(1/popt[3]))
print(" s0 =  "+str(popt[0]))
print(" a  = "+str(popt[1]*popt[3]))
print(" b  =  "+str(popt[2]*popt[3]))

#old version: s0=1.99 K=1.162 a=1.012 b=-0.258

# %% q exponent law plot

plt.figure(figsize=(7,4))
plt.errorbar(bin_sizes, exps, yerr=3*errs, color='firebrick',fmt='o', label='$q_{  \ j}\ $ slope fit',markerfacecolor='white',linewidth=1.5,markeredgewidth=1.5)
plt.xscale('log')
plt.plot(bin_sizes,dust_radial_exponent(bin_sizes,*popt),linewidth=2.,linestyle='--',color='mediumseagreen', label='$q (s)\ $ empyrical law')
#plt.plot(bin_sizes,dust_radial_exponent_p(bin_sizes),linewidth=3., label='$q\ (s)$ fit')

plt.xlabel("dust size $s$ [cm]")
plt.ylabel(r"radial distribution exponent  $q=\frac{d \ \log N_P}{d \ \log r}$")
plt.legend(loc='upper center',fancybox=False,framealpha=1.)
plt.grid(which='both')
plt.savefig('output_plots/q_fit.pdf',bbox_inches='tight')


# %% particle counts

cmap=nice_plots()

plt.figure(figsize=(10,5))
for i in range(N_bins):

    j= N_bins -i -1
    formatted_size = f"{bin_sizes[j]:.3f}"
    plt.plot(time[::100],counts[::100,j],label=formatted_size,c=cmap(i/N_bins),alpha=1)
plt.xlabel('time $t$ [$2\pi / \Omega_p$]')
plt.ylabel('Nuber of particles per size bin $N_P$')
plt.grid()
#
plt.legend(fontsize=12,title='particle size $s$ [cm]',ncols=5, bbox_to_anchor=(0.5, -0.2), loc='upper center', fancybox=False, framealpha=0.)

plt.savefig('output_plots/particle_counts.pdf',bbox_inches='tight')

# %% Stokes Number

output_dir = r'./problem/out/'

cmap=nice_plots()

def get_stokes(step):

    
    AU = 1.496e11 # Astronomical unit in meters
    MSUN = 1.98e30 # Solar mass in kg
    YEAR = 3.1557e7 # Year in seconds

    rp=1 * AU
    M_star= 1 * MSUN
    t_sim = 1 * YEAR

    delta_r= 0.01 * AU
    
    G = 6.6743e-11 # in kms units

    x, y, size_part, size_bins, time, x_planet, y_planet, t_planet, tau_part = positions (step)
    r=np.sqrt(x**2+y**2) * AU
    tau_part = tau_part * YEAR

    OmegaK=np.sqrt(G*M_star/(rp**3))

    stokes_part = tau_part*OmegaK

    r_filter = (r>rp-delta_r) & (r<rp+delta_r)

    stokes_avg = np.zeros(len(size_bins))
    stokes_err = np.zeros(len(size_bins))

    for i in range(len(size_bins)):

        size_filter = (size_part == size_bins[i])
        stokes_filtered = stokes_part[r_filter & size_filter]
        stokes_avg[i] = np.average(stokes_filtered)
        stokes_err = 0.5*(max(stokes_filtered)-min(stokes_filtered))

    size_bins = size_bins * AU *100

    return size_bins, stokes_avg, stokes_err

size, stokes, stokes_err = get_stokes(2)

plt.figure()
plt.plot(size,stokes)
plt.grid(which='both')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Particle size [cm]')
plt.ylabel(r'Stokes number $\mathcal{S}$')

def plot_stokes(step,cmap):

    
    AU = 1.496e11 # Astronomical unit in meters
    MSUN = 1.98e30 # Solar mass in kg
    YEAR = 3.1557e7 # Year in seconds

    rp=1 * AU
    M_star= 1 * MSUN
    t_sim = 1 * YEAR

    delta_r= 1 * AU
    
    G = 6.6743e-11 # in kms units

    x, y, size_part, size_bins, time, x_planet, y_planet, t_planet, tau_part = positions (step)
    r=np.sqrt(x**2+y**2) * AU
    tau_part = tau_part * YEAR

    OmegaK=np.sqrt(G*M_star/(rp**3))
    Omega_part = np.sqrt(G*M_star/(r**3))
    
    stokes_part = tau_part*Omega_part

    r_filter = (r>rp-delta_r) & (r<rp+delta_r)

    for i in range(len(size_bins)):

        size_filter = (size_part == size_bins[i])
        stokes_filtered = stokes_part[r_filter & size_filter]
        radius_filtered = r[r_filter & size_filter]/AU
        plt.scatter(stokes_filtered, radius_filtered, s=1,label=f'size {size_bins[i]:.3f} cm', alpha=0.5, c=cmap(i/len(size_bins)))

plt.figure()
plot_stokes(2,cmap)
plt.xscale('log')
plt.grid(which='both')
plt.xlabel(r'Stokes number $\mathcal{S}$')
plt.ylabel('Radial position [AU]')
# %%
