#%%
import numpy as np
import pyPLUTO.pload as pp 
import pyPLUTO.Image as img
import particles_read as pr
import matplotlib.pyplot as plt
import astropy.units as u
from scipy.stats import norm
from scipy.optimize import curve_fit
from pathlib import Path
from matplotlib import font_manager as fm, cycler
import os

global output_dir

output_dir = r'.problem/out/'

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

    # --- Math font settings ---
    if use_math:
        plt.rcParams.update({
            "mathtext.fontset": "custom",
            "mathtext.rm": "EB Garamond",
            "mathtext.it": "EB Garamond:italic",
            "mathtext.bf": "EB Garamond:bold",
        })

    print("✅ EB Garamond style configured for Matplotlib.")

def positions(step):

    # import dust data
    P = pr.ReadPartData(step)

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

    return x_part, y_part, size_part, size_bins , time, x_planet, y_planet, t_planet

def dust_plot(step,S):

    x, y, size_part, size_bins, time, x_planet, y_planet, t_planet = positions(step)
    
    r = np.sqrt(y**2+x**2)
    angle = np.arctan2(y,x)

    r_planet = np.sqrt(x_planet**2+y_planet**2)
    angle_planet = np.arctan2(y_planet,x_planet)
    # plt.scatter(angle, radius, s=0.2, c=size_part, cmap='viridis')
    # plt.scatter(angle_planet, r_planet, s=10, c='black')

    angle_relative = angle -angle_planet
    angle_relative[angle_relative>np.pi] -= 2*np.pi
    angle_relative[angle_relative<-np.pi] += 2*np.pi

    plt.figure()
    size_index = (size_part==size_bins[S])
    plt.scatter(angle_relative[size_index], r[size_index], s=0.2, c=size_part[size_index], cmap='brg')
    plt.scatter(0, r_planet, s=10, c='black')
    plt.xlim(-np.pi,np.pi)

    plt.figure()
    plt.hist(r[size_part==size_part[S]], bins=100)

def dust_hist(step,S):

    plt.figure()
    for i in range(len(step)):
        x, y, size_part, size_bins, time, x_planet, y_planet, t_planet = positions(step[i])
        r = np.sqrt(y**2+x**2)
        #r_planet = np.sqrt(x_planet**2+y_planet**2)

        plt.hist(r[size_part==size_part[S]], alpha=0.8, bins=200,label='time: '+str(time))

    plt.legend()
    plt.show()

def smoothing_length(r_planet):
    
    smoothing_factor=1.
    m_ratio = 1e-4 #planet to star mass ratio

    hill_radius = r_planet * (m_ratio/3)**(1/3)
    d_smooth = smoothing_factor*hill_radius

    return d_smooth

def calculate_x_force(step):

    x_part, y_part, size_part, time, r_planet, angle_planet, t_planet =dust_positions(step)

    x_planet = r_planet * np.sin(angle_planet)
    y_planet = r_planet * np.cos(angle_planet)

    x_part_rel = x_part - x_planet
    y_part_rel = y_part - y_planet

    d_smooth2 = smoothing_length(r_planet)**2

    d2 = x_part_rel**2 +y_part_rel**2 + d_smooth2
    force_x_part = x_part_rel / d2**1.5

    sizes = np.sort(np.unique(size_part))

    force_size = np.zeros(len(sizes))

    #for each size bin sum the force

    for i, size in enumerate(sizes):
        size_mask = size_part == size
        force_size[i] = np.sum(force_x_part[size_mask])
        
    return sizes, force_size, time, r_planet

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

def mov_stats(X, window_size):

    half = window_size // 2
    X_avg = np.zeros_like(X, dtype=float)
    X_err = np.zeros_like(X, dtype=float)
    
    for i in range(len(X)):
        start = max(0, i - half)
        end = min(len(X), i + half + 1)
        window = X[start:end]
        X_avg[i] = np.mean(window, axis=0)
        X_err[i] = np.std(window, axis=0)/np.sqrt(len(window))
    
    return X_avg, X_err

def plot_dist_gaussian(data):

    plt.figure()
    counts, bins, _ = plt.hist(data, bins='auto', density=True, alpha=0.6, color='mediumturquoise')

    # Fit a normal distribution to the data
    avg, sigma = norm.fit(data)
    # Generate the Gaussian curve
    x = np.linspace(bins[0], bins[-1], 100)
    pdf = norm.pdf(x, avg, sigma)

    err = sigma/np.sqrt(len(data))

    # Plot the fitted Gaussian
    plt.plot(x, pdf, 'darkred', linewidth=2, label=f'Fit: μ={avg:.2e}, σ={sigma:.2e}',alpha=0.5)
    plt.plot([avg, avg], [0, norm.pdf(avg, avg, sigma)], color= 'teal',linestyle='--')
    plt.legend()
    plt.xlabel("Value")
    plt.ylabel("Distribution density")
    plt.title("Force distribution")
    plt.show("")

    return avg, sigma, err

def radial_distribution(S,step_start,step_end, plotting=True):

    from scipy.optimize import curve_fit

    N_steps = step_end - step_start
    radii=np.zeros(N_steps)

    for i in range(N_steps):

        step= step_start +i
        x, y, size_part, size_bins, time, x_planet, y_planet, t_planey = positions (step)

        size_index = (size_part==size_bins[S])
        r= np.sqrt(x[size_index]**2 +y[size_index]**2)
        radii=np.append(radii,r)

    # fit power law

    rmax = 1.8*0.9
    rmin = 0.5*1.1
    r_outer = rmax-0.1*(rmax-rmin)

    if plotting:
        plt.figure()
        counts, bins, _ = plt.hist(radii[(radii<r_outer) & (radii>rmin)], bins=100, color='teal', density=True)
    else:
        counts, bins = np.histogram(radii[(radii<r_outer) & (radii>rmin)], bins=100, density=True)

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

#%%

dust_hist(step=[0,9],S=2)

#%%
size=1
t_stat = 0

nice_plots()

time, force, counts, N_bins, bin_sizes = read_force()

#dust slope calculation

exps = np.zeros(N_bins)
errs = np.zeros(N_bins)

for i in range(N_bins-1):
    exps[i], errs[i] = radial_distribution(i,1,9,plotting=True)

def dust_radial_exponent(x,x0,a,b,n):
    #for model fitting
    
    u=x/x0
    s1=u**(a*n)
    s2=u**(b*n)

    exponent= np.log10(s1+s2)/n

    return exponent

def dust_radial_exponent_p(x):
    #fitted model

    u=x/1.99
    return 1.162*np.log10(u**1.012+u**-0.258)

popt, pcov = curve_fit(dust_radial_exponent,bin_sizes,exps,[2,-0.3,1.2,0.86],errs)

print(popt)
print(np.linalg.cond(pcov))

#%%
plt.figure()
plt.errorbar(bin_sizes, exps, yerr=errs, color='darkgreen',fmt='o', alpha=0.5,label='simulation ($\pm$ 1 std)')
plt.plot(bin_sizes,dust_radial_exponent_p(bin_sizes),linewidth=3., label='$q\ (s)$ fit')

plt.xscale('log')
plt.xlabel("dust size $s$ [cm]")
plt.ylabel(r"radial distribution exponent  $\frac{d \ \log N_p}{d \ \log r}$")
plt.legend(loc='upper center',fancybox=False,framealpha=1.)
plt.title('Steady state dust exponent')
plt.grid(which='both')

# %% torque calculation

#time, force, counts, N_bins, bins = read_force()

# avg_force = np.mean(force[time>t_stat,size])

# avg_counts = np.mean(counts[time>t_stat,:],axis=0)

# avg_force, sigma_force, err = plot_dist_gaussian(force[time>t_stat,size])

# plt.show()
# print(f"Average force: {avg_force:.3e} +/- {3*err:.3e}")
# print(f"Force std: {sigma_force:.3e}")

# %% particle counts

plt.figure()
for i in range(N_bins):
    plt.plot(time,counts[:,i],label=str(bin_sizes[i]),c=plt.cm.viridis(i/N_bins))
plt.xlabel('Time')
plt.legend(loc='center left',bbox_to_anchor=(1, 0.5))
plt.ylabel('Nuber of particles in domain')

plt.figure()
plt.semilogx(bin_sizes,avg_counts,'o-')

# %%
