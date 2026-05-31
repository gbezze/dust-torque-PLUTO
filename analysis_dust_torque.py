#%% Setup
import numpy as np
import pyPLUTO.pload as pp 
import pyPLUTO.Image as img
#import particles_read as pr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import astropy.constants as const
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.interpolate import interp2d, interp1d
from scipy.integrate import quad
from pathlib import Path
from matplotlib import font_manager as fm, cycler
import os
import warnings
import mpl_scatter_density
import matplotlib.patches as patches
import pickle as pk
import pandas as pd



global output_dir

output_dir = r'./problem/out/'

warnings.filterwarnings("ignore", category=DeprecationWarning)

def nice_plots():
    
    font_dir = r"/home/bezze/ebgaramond"
    palette = None
    font_size = 14
    use_math = True

    font_dir = Path(font_dir)
    if not font_dir.exists():
        raise FileNotFoundError(f"Font directory not found: {font_dir}")

    # Register fonts
    for font_file in font_dir.glob("EBGaramond*.ttf"):
        fm.fontManager.addfont(str(font_file))

    # Set base font
    plt.rcParams["font.family"] = "EB Garamond"
    plt.rcParams["font.size"] = font_size
    plt.rcParams["axes.titleweight"] = "bold"

    # Color palette
    if palette is None:
        palette = ["firebrick", "mediumseagreen", "darkorchid", "gold", "deepskyblue"]
    plt.rcParams["axes.prop_cycle"] = cycler(color=palette)

    # Grid color
    plt.rcParams['grid.color'] = 'black'
    plt.rcParams['grid.alpha'] = 0.2
    plt.rcParams['grid.linewidth'] = 0.7

    #legend
    plt.rcParams['legend.fancybox'] = False

    # Math font settings
    if use_math:
        plt.rcParams.update({
            "mathtext.fontset": "custom",
            "mathtext.rm": "EB Garamond",
            "mathtext.it": "EB Garamond:italic",
            "mathtext.bf": "EB Garamond:bold",
        })

    print("EB Garamond configured")

def ReadPartData(ns):

    def NStepStr(ns):
        nsstr = str(ns)
        while len(nsstr) < 4:
            nsstr = '0'+nsstr
        return nsstr
    
    nstepstr = NStepStr(ns)
    fname = output_dir+"particles."+nstepstr+".partdbl"
    h_lines = 0 
    val_dict = {}

    #READ HEADER. 
    with open(fname, "rb") as f:
        for line in f:
            if h_lines < 13:
                if h_lines > 0 and h_lines < 13:
                    val_dict.update({line.split()[1].decode('utf8'):[i.decode('utf8') for i in line.split()[2:]]})
                h_lines += 1

    fdata = open(fname, "rb")
 
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

def get_sym_data(printing = False):

    global M_ratio, sigma0, rid, rod, rb, timestep

    planet_data_path = os.path.join(output_dir, 'nbody.out')
    logfile_path = os.path.join(output_dir, 'pluto.0.log')

    with open(planet_data_path,'r') as f:

        lines = f.readlines()
        M_ratio = np.array([float(lines[1].split()[2].strip())])

    with open(logfile_path,'r') as f:

        lines = f.readlines()
        for line in lines:
            if 'SIGMA_0' in line:
                sigma0 = float(line.split()[2].strip())
            if 'DUST_INJECTION_RATIO' in line:
                rb_ratio = float(line.split()[2].strip())
            else:
                rb_ratio = 0.1
            if 'DAMPING_INNER' in line:
                damping_inner = float(line.split()[2].strip())   
            if 'DAMPING_OUTER' in line:
                damping_outer = float(line.split()[2].strip())
            if 'X1-grid' in line:
                ri = float(line.split()[3].strip())
                ro = float(line.split()[6].strip())

        lines = lines[-20:-1]

        for line in lines:
            if r'dt =' in line:
                timestep = float(line.split()[6].strip()[:-1])/2

    rid = ri * damping_inner
    rod = ro * damping_outer
    rb = rod - (rod-rid)*rb_ratio

    if(printing):
        print("")
        print("Simulation parameters:")
        print(f"Planet-to-star mass ratio: {M_ratio[0]:.2e}")
        print(f"Gas surface density at 1 AU: {sigma0:.2e} g/cm2")
        print(f"Dust inner radius: {rid:.2f} AU")
        print(f"Dust outer radius: {rod:.2f} AU")
        print(f"Dust injection boundary: {rb:.2f} AU")
        print(f"Timestep: {timestep:.2e} OmegaK^-1")
        print("")

def positions(step):

    # import dust data
    P = ReadPartData(step)

    r_part = np.array(P['x1']).flatten()
    phi_part = np.array(P['x2']).flatten()
    tau_part = np.array(P['tau_s']).flatten()
    size_part = np.array(P['radius']).flatten()

    t_part = float(P['time'][0])

    #import planet position

    planet_data_path = os.path.join(output_dir, 'nbody_coordinates.dat')

    with open(planet_data_path,'r') as f:

        lines = f.readlines()
        times = np.array([float(line.split()[1]) for line in lines[9:]])
        time_index = np.argmin(np.abs(times - t_part))

        x_planet = float(lines[time_index].split()[2])
        y_planet = float(lines[time_index].split()[3])
        t_planet = float(lines[time_index].split()[1])

    x_part = r_part * np.cos(phi_part) 
    y_part = r_part * np.sin(phi_part)

    truncate = lambda x, n: np.where(x==0, 0, np.trunc(x/10**(np.floor(np.log10(np.abs(x)))-n+1))*10**(np.floor(np.log10(np.abs(x)))-n+1))
    
    size_part = truncate(size_part,7)
    size_bins = np.sort(np.unique(size_part))

    r_planet = np.sqrt(x_planet**2 + y_planet**2)
    angle_planet = np.arctan2(y_planet,x_planet)

    time_difference = t_part - t_planet

    delta_angle = 2*np.pi * (time_difference + timestep) * np.sqrt(1/(r_planet**3))

    x_planet = r_planet * np.cos(angle_planet + delta_angle)
    y_planet = r_planet * np.sin(angle_planet + delta_angle)

    # calculate individual particle tangential force on planet

    return x_part, y_part, size_part, size_bins , t_part, x_planet, y_planet, t_planet, P

def dust_plot(step, S, zoomangle=1., zoomr=1.):

    x, y, size_part, size_bins, time, x_planet, y_planet, t_planet, P = positions(step)
    
    r = np.sqrt(y**2+x**2)
    angle = np.arctan2(y,x)

    r_planet = np.sqrt(x_planet**2+y_planet**2)
    angle_planet = np.arctan2(y_planet,x_planet)

    angle_relative = angle -angle_planet
    angle_relative[angle_relative>np.pi] -= 2*np.pi
    angle_relative[angle_relative<-np.pi] += 2*np.pi

    deltar = 0.5/zoomr #normalized to r_p
    deltaangle = 1/zoomangle #normalzied to pi

    size_index = np.zeros(len(r),dtype=bool)

    try:
        for i in range(len(S)):
            size_index[size_part==size_bins[S[i]]] = True
    except:
        size_index = (size_part==size_bins[S])

    # plt.scatter(angle_relative[size_index], r[size_index], s=0.2, c=size_part[size_index], cmap='brg',alpha=0.7)
    plt.scatter(angle_relative[size_index], r[size_index], s=0.5, c='darkred',alpha=1)
    plt.scatter(0, r_planet, s=10, c='black')
    plt.xlim(-deltaangle*np.pi,deltaangle*np.pi)
    plt.ylim(r_planet*(1-deltar),r_planet*(1+deltar))
    plt.xlabel(r'$\phi - \phi_p$ [rad]')
    plt.ylabel(r'$r$ [AU]')

    r_h = 0.1 * r_planet * (M_ratio/3)**(1/3)

    theta = np.linspace(0, 2*np.pi, 400)
    a_hill = (r_h/r_planet) * np.cos(theta)
    r_hill = r_planet + r_h * np.sin(theta)

    #plt.plot(a_hill, r_hill, color='forestgreen', linewidth=1.5)
    # plt.text(0, r_planet + r_h, r'$r_H$', color='forestgreen', verticalalignment='bottom', horizontalalignment='left')

    # plt.figure()
    # plt.hist(r[(size_part==size_part[S])], bins=100)
    # plt.xlabel(r'$r$ [AU]')
    # plt.ylabel('Particles count')
    # plt.figure()
    #plt.hist(angle_relative[(size_part==size_part[S])&(np.abs(angle_relative)<deltaangle*np.pi)&(np.abs(r-r_planet)<deltar*r_planet)], bins=100)

    #return P

def dust_densplot(step, S, zoomangle=1., zoomr=1.,maxcounts=4):

    x, y, size_part, size_bins, time, x_planet, y_planet, t_planet, P = positions(step)

    r = np.sqrt(y**2+x**2)
    angle = np.arctan2(y,x)

    r_planet = np.sqrt(x_planet**2+y_planet**2)
    angle_planet = np.arctan2(y_planet,x_planet)

    angle_relative = angle -angle_planet
    angle_relative[angle_relative>np.pi] -= 2*np.pi
    angle_relative[angle_relative<-np.pi] += 2*np.pi

    deltar = 0.5/zoomr #normalized to r_p
    deltaangle = 1/zoomangle #normalzied to pi

    size_index = np.zeros(len(r),dtype=bool)

    try:
        for i in range(len(S)):
            size_index[size_part==size_bins[S[i]]] = True
    except:
        size_index = (size_part==size_bins[S])

    cmap = plt.get_cmap('magma_r')
    colors = cmap(np.linspace(0, 1, 100))
    colors[0] = [1, 1, 1, 1]  # white
    cmap_white=mcolors.LinearSegmentedColormap.from_list(
        f"custom_white", colors
    )
    

    fig = plt.figure(figsize=(5,4))
    ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
    norm = mcolors.Normalize(vmin=0, vmax=maxcounts)
    dust_density=ax.scatter_density(angle_relative[size_index], r[size_index], cmap=cmap_white, norm=norm)
    plt.scatter(0, r_planet, s=10, c='black')
    plt.xlim(-deltaangle*np.pi,deltaangle*np.pi)
    plt.ylim(r_planet*(1-deltar),r_planet*(1+deltar))
    # plt.xlabel(r'$\phi - \phi_p$ [rad]')
    # plt.ylabel(r'$r$ [AU]')
    
    plt.ylabel('distance from star [AU]')
    plt.xlabel('angle relative to planet [rad]')


def dust_hist(step,S):

    plt.figure()
    for i in range(len(step)):
        x, y, size_part, size_bins, time, x_planet, y_planet, t_planet = positions(step[i])
        r = np.sqrt(y**2+x**2)
        #r_planet = np.sqrt(x_planet**2+y_planet**2)

        plt.hist(r[size_part==size_part[S]], alpha=0.8, bins=200,label='time: '+str(time))

    plt.legend()
    plt.show()
    
    smoothing_factor=1.
    m_ratio = 1e-4 #planet to star mass ratio

    hill_radius = r_planet * (m_ratio/3)**(1/3)
    d_smooth = smoothing_factor*hill_radius

    return d_smooth

def read_force():
    
    force_data_path = os.path.join(output_dir, 'dust_force.dat')

    with open(force_data_path,'r') as f:

        lines = f.readlines()
        N_bins = int(lines[0].split(": ")[1].strip())
        N_times = len(lines)-6

        time=np.zeros(N_times)
        force=np.zeros((N_times,N_bins))
        force=np.zeros((N_times,N_bins))
        counts=np.zeros((N_times,N_bins))

        bins=np.zeros(N_bins)
        bins = np.array(lines[2].split()[3:N_bins+3],dtype=float)

        for i in range(N_times):
            time[i]=lines[i+6].split()[0]
            force[i,:]=lines[i+6].split()[1:N_bins+1]
            counts[i,:]=lines[i+6].split()[N_bins+1:2*N_bins+1]

    return time, force, counts, N_bins, bins

def fit_gaussian(data):

    def gaussian(x, mu, sigma):
        A = (2 * np.pi * sigma**2)**-0.5
        return A * np.exp(-0.5 * ((x - mu) / sigma)**2)

    # Fit a normal distribution to the data
    #use norm.fit first (fast and robust)
    avg, sigma = norm.fit(data)

    #use curve_fit (estimates covariance)
    counts, bins = np.histogram(data, bins='auto', density=True)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    popt, pcov= curve_fit(gaussian, bin_centers, counts, [avg,sigma], maxfev=10000)

    avg=popt[0]
    sigma=popt[1]

    # Generate the Gaussian curve
    x = np.linspace(bins[0], bins[-1], 100)
    pdf = norm.pdf(x, avg, sigma)

    #err = sigma/np.sqrt(len(data))
    err=3*np.sqrt(pcov[1,1])

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

def dust_radial_exponent_p(size):
    #fitted model, size in cm

    u=size/2.827
    return 2.246*np.log10(u**0.627+u**-0.066)

def torque_0():

    #all in kg m s

    # planet orbit
    G = 6.674e-11
    rp = 1.496e11
    M_star = 1.989e30
    Omega_k = np.sqrt(G*M_star/(rp**3))

    # disk properties
    h=0.05
    Sigma_0 = sigma0 * 10 #[kg/m2]

    torque= M_ratio**2 * rp**4 * Omega_k**2 * Sigma_0 / h**2

    return torque

def torque_factor_density(size):
    #computes the torque factor to have a prescribed local dust/gas ratio

    G = 6.6743e-11 # [m3 / kg s2]
    AUtom = 1.496e11 # Astronomical unit in meters (AU * AUtom = m)

    q=dust_radial_exponent_p(size)

    rp_au = 1 * AUtom
    rid_au = rid * AUtom
    rod_au = rod * AUtom
    rb_au = rb * AUtom

    Sigma_gas = sigma0*10 #[kg / m2]
    
    dust_gas_ratio = 0.01

    M_dust = 2 * np.pi * Sigma_gas * dust_gas_ratio * (rb_au**(q+2) - rid_au**(q+2)) / ((q+2)*rp_au**q)

    # PLUTO calculates y/r^3 in units of 1/AU^2 so we need to convert it
    conv_factor=AUtom**2

    torque_factor=G*M_dust*rp_au*conv_factor
    
    return torque_factor

def torque_gas(n):
    
    torque = -(1.364-0.541*n)*torque_0()

    return torque

def autocorrelation(f):

    n = len(f)    
    f -= np.mean(f)
    
    # FFT-based autocorrelation
    F = np.fft.fft(f, n=2*n)
    acf = np.fft.ifft(F * np.conjugate(F))[:n].real

    # Unbiased normalization
    acf /= (np.arange(n, 0, -1))
    
    # Normalize so that acf[0] = 1
    acf /= acf[0]

    # Sokal's method to estimate integrated autocorrelation time
    c = 5
    tau = 0.5  # initial value

    for W in range(1, n):
        # compute τint(W) = 1/2 + sum_{t=1}^W ρ(t)
        tau_new = 0.5 + np.sum(acf[1:W])

        # stopping criterion: W > c * τint(W)
        if W > c * tau_new:
            tau = tau_new
            break
        tau = tau_new

    return acf, tau

def stats(data):

    # compute average, std and error considering autocorrelation time

    avg = np.mean(data)
    sigma = np.std(data)

    n = len(data)
    _, n_acf = autocorrelation(data)
    N_effective = n / (2*n_acf +1)

    err = 3*sigma/np.sqrt(N_effective)

    return avg, sigma, err

def torque_BLP18(stokes):

    #data from the table
    M_ratio_values = np.array([0.333, 0.486, 0.709, 1.03, 1.51, 2.2, 3.22, 4.69, 6.85, 10.0])
    stokes_values = np.array([0.010, 0.014, 0.021, 0.030, 0.043, 0.062, 0.089, 0.127, 0.183, 0.264, 0.379, 0.546, 0.785, 1.129])
    torque_values = np.array([
        [0.283, 0.309, 0.264, 0.201, 0.142, 0.093, 0.056, 0.029, 0.006, -0.010],
        [0.371, 0.394, 0.415, 0.333, 0.239, 0.160, 0.098, 0.053, 0.017, -0.008],
        [0.739, 0.720, 0.618, 0.475, 0.382, 0.276, 0.172, 0.095, 0.041, -0.002],
        [1.638, 1.332, 1.102, 0.851, 0.625, 0.434, 0.275, 0.185, 0.087, 0.021],
        [3.254, 2.614, 1.961, 1.476, 1.078, 0.757, 0.503, 0.316, 0.174, 0.085],
        [5.776, 4.579, 3.428, 2.461, 1.790, 1.266, 0.850, 0.546, 0.326, 0.167],
        [9.672, 7.297, 5.313, 3.760, 2.620, 1.872, 1.301, 0.880, 0.562, 0.304],
        [13.684, 10.191, 7.337, 5.145, 3.471, 2.374, 1.725, 1.209, 0.805, 0.483],
        [16.058, 11.276, 7.866, 5.231, 3.178, 1.372, 1.740, 1.493, 1.073, 0.701],
        [-1.521, -18.135, -8.446, -0.425, 1.971, 1.746, 2.237, 1.823, 1.351, 0.910],
        [-32.124, -7.348, 2.307, 4.517, 4.752, 3.232, 2.803, 2.108, 1.499, 1.016],
        [5.523, 8.831, 9.054, 7.882, 5.860, 3.741, 2.716, 1.953, 1.398, 0.986],
        [10.492, 9.803, 7.187, 5.859, 4.314, 3.278, 2.465, 1.821, 1.326, 0.943],
        [11.861, 9.572, 7.850, 5.807, 4.538, 3.367, 2.463, 1.799, 1.299, 0.916]])

    #interpolate in log-log space
    log_M_ratio = np.log10(M_ratio_values)
    log_stokes = np.log10(stokes_values)
    interp_func = interp2d(log_M_ratio, log_stokes, torque_values, kind='linear')
    torque = interp_func(np.log10(M_ratio), np.log10(stokes))

    return torque

def torque_rescale_factor(size, dust_gas_ratio):
    
    G = 6.6743e-11 # [m3 / kg s2]
    AU = 1.496e11 # Astronomical unit in meters
    M_star = 1.989e30 # Solar mass in kg
    Sigma_gas = sigma0 *10 #[kg / m2]

    M_planet = M_ratio * M_star

    q=dust_radial_exponent_p(size)

    rp_m = 1 * AU
    rid_m = rid * AU
    rod_m = rod * AU
    rb_m = rb * AU

    M_dust = 2 * np.pi * Sigma_gas * dust_gas_ratio * (rb_m**(q+1) - rid_m**(q+1)) / ((q+1)* (rp_m**(q-1)))

    torque_rescale_factor = G * M_dust * M_planet /rp_m

    return torque_rescale_factor

def get_stokes(step=1):

    AU = 1.496e11 # Astronomical unit in meters
    MSUN = 1.98e30 # Solar mass in kg
    YEAR = 3.1557e7 # Year in seconds

    rp=1 * AU
    M_star= 1 * MSUN
    t_sim = 1 * YEAR 

    delta_r= 2*1.01 * rp * np.sqrt(M_ratio/0.05) # corotation zone
    
    G = 6.6743e-11 # in kms units



    x, y, size_part, size_bins, time, x_planet, y_planet, t_planet, tau_part = positions (step)
    r=np.sqrt(x**2+y**2) * AU
    tau_part = np.array(tau_part['tau_s']) * YEAR  # Replace 'tau_s' with the correct key if needed

    OmegaK=np.sqrt(G*M_star/(rp**3))

    stokes_part = tau_part*OmegaK

    r_filter = (r>rp-delta_r) & (r<rp+delta_r)

    stokes_avg = np.zeros(len(size_bins))
    stokes_err = np.zeros(len(size_bins))

    for i in range(len(size_bins)):

        size_filter = (size_part == size_bins[i])
        stokes_filtered = stokes_part[r_filter & size_filter]
        stokes_avg[i] = np.average(stokes_filtered)
        stokes_err[i] = 0.5*(max(stokes_filtered)-min(stokes_filtered))
        #stokes_err[i] = 3*np.std(stokes_filtered-stokes_avg[i])

    size_bins = size_bins * AU *100

    return size_bins, stokes_avg, stokes_err

def plot_BLP18(size1=0.2,size2=21):

    #data from the table
    M_ratio_values = np.array([0.333, 0.486, 0.709, 1.03, 1.51, 2.2, 3.22, 4.69, 6.85, 10.0])*3e-6
    stokes_values = np.array([0.010, 0.014, 0.021, 0.030, 0.043, 0.062, 0.089, 0.127, 0.183, 0.264, 0.379, 0.546, 0.785, 1.129])
    torque_values = np.array([
        [0.283, 0.309, 0.264, 0.201, 0.142, 0.093, 0.056, 0.029, 0.006, -0.010],
        [0.371, 0.394, 0.415, 0.333, 0.239, 0.160, 0.098, 0.053, 0.017, -0.008],
        [0.739, 0.720, 0.618, 0.475, 0.382, 0.276, 0.172, 0.095, 0.041, -0.002],
        [1.638, 1.332, 1.102, 0.851, 0.625, 0.434, 0.275, 0.185, 0.087, 0.021],
        [3.254, 2.614, 1.961, 1.476, 1.078, 0.757, 0.503, 0.316, 0.174, 0.085],
        [5.776, 4.579, 3.428, 2.461, 1.790, 1.266, 0.850, 0.546, 0.326, 0.167],
        [9.672, 7.297, 5.313, 3.760, 2.620, 1.872, 1.301, 0.880, 0.562, 0.304],
        [13.684, 10.191, 7.337, 5.145, 3.471, 2.374, 1.725, 1.209, 0.805, 0.483],
        [16.058, 11.276, 7.866, 5.231, 3.178, 1.372, 1.740, 1.493, 1.073, 0.701],
        [-1.521, -18.135, -8.446, -0.425, 1.971, 1.746, 2.237, 1.823, 1.351, 0.910],
        [-32.124, -7.348, 2.307, 4.517, 4.752, 3.232, 2.803, 2.108, 1.499, 1.016],
        [5.523, 8.831, 9.054, 7.882, 5.860, 3.741, 2.716, 1.953, 1.398, 0.986],
        [10.492, 9.803, 7.187, 5.859, 4.314, 3.278, 2.465, 1.821, 1.326, 0.943],
        [11.861, 9.572, 7.850, 5.807, 4.538, 3.367, 2.463, 1.799, 1.299, 0.916]])

    #interpolate in log-log space
    log_M_ratios = np.log10(M_ratio_values)
    log_stokes = np.log10(stokes_values)

    size_bins, stokes_avg, stokes_err = get_stokes()
    interp_size=interp1d(size_bins, stokes_avg, kind='linear', fill_value='extrapolate')
    interp_stokes=interp1d(stokes_avg, size_bins, kind='linear', fill_value='extrapolate')
    interp_stokes_plus=interp1d(stokes_avg+stokes_err, size_bins, kind='linear', fill_value='extrapolate')
    interp_stokes_minus=interp1d(stokes_avg-stokes_err, size_bins, kind='linear', fill_value='extrapolate')
    stokes1=interp_size(size1)
    stokes2=interp_size(size2)

    stokes_included = stokes_values[(stokes_values > stokes1) & (stokes_values <stokes2)]
    
    interp_torques = interp2d(log_M_ratios, log_stokes, torque_values, kind='linear')
    torques = interp_torques(np.log10(M_ratio), np.log10(stokes_included))
    sizes_effective = interp_stokes(stokes_included)
    sizes_effective_plus = interp_stokes_plus(stokes_included)
    sizes_effective_minus = interp_stokes_minus(stokes_included)

    print(f'Showing BLP18 torque for Mp/Mstar = {M_ratio[0]:.2e}')
    
    plt.plot(sizes_effective, torques.flatten(), '--', label='BLP18', color='forestgreen',alpha=1)
    plt.fill_betweenx(torques.flatten(), sizes_effective_minus, sizes_effective_plus, color='forestgreen', alpha=0.3, linewidth = 0)

def plot_sym_torque(Color='steelblue',Label='simulation'):

    get_sym_data(printing=True)

    time, torque_normalized, counts, N_bins, bins = read_force()

    print(time[-1])

    tstat = 200

    dust_gas_ratio = 0.01

    torque_d=np.zeros(N_bins)
    torque_d_err=np.zeros(N_bins)

    for i in range(N_bins):

        torque_factor=torque_rescale_factor(bins[i],dust_gas_ratio)
        xi_avg, xi_sigma, xi_err = stats(torque_normalized[time>tstat,i]*torque_factor)
        torque_d[i]= xi_avg
        torque_d_err[i] = xi_err

    plt.plot(bins, torque_d/torque_0(), '.-', label=Label, color=Color,markerfacecolor=Color,linewidth=1.5,markeredgewidth=1.5)
    plt.fill_between(bins, (torque_d - torque_d_err)/torque_0(), (torque_d + torque_d_err)/torque_0(), color=Color, alpha=0.2, linewidth=0)

def stokes_ticks(ax):

    sizes, stokes, stokes_err = get_stokes(1)

    def size_to_stokes(size):
        return np.interp(size, sizes, stokes)

    def stokes_to_size(stoke):
        return np.interp(stoke, stokes, sizes)

    secax = ax.secondary_xaxis(
        'top',
        functions=(size_to_stokes, stokes_to_size)
    )

    secax.tick_params(axis='x', pad=0)

    def hide_first_n_labels(ax, n=1, axis='x', color='white', pad=0.001):
        """
        Hide the first n tick labels on a given axis by drawing white boxes over them.

        Parameters:
            ax      : matplotlib.axes.Axes object (primary or secondary)
            n       : number of first labels to hide
            axis    : 'x' or 'y'
            color   : color of the box (usually same as background)
            pad     : padding around the label (optional, in figure coordinates)
        """
        # Force draw to get the label positions
        ax.figure.canvas.draw()
        if axis.lower() == 'x':
            labels = ax.get_xticklabels()
        elif axis.lower() == 'y':
            labels = ax.get_yticklabels()
        else:
            raise ValueError("axis must be 'x' or 'y'")

        # Draw white rectangles over the first n labels
        renderer = ax.figure.canvas.get_renderer()
        for label in labels[:n]:
            bbox = label.get_window_extent(renderer=renderer)
            bbox_fig = bbox.transformed(ax.figure.transFigure.inverted())
            rect = patches.FancyBboxPatch(
                (bbox_fig.x0 - pad, bbox_fig.y0 - pad),
                bbox_fig.width + 2*pad,
                bbox_fig.height + 2*pad,
                boxstyle="square,pad=0",
                transform=ax.figure.transFigure,
                color=color,
                zorder=5
            )
            ax.figure.patches.append(rect)

    hide_first_n_labels(secax)

    secax.set_xlabel(r'particle Stokes number $\mathcal{S}$',labelpad=5)

nice_plots()

# %% dust scatter plot
output_dir = r'./problem/out_D3/'

get_sym_data()
filenum=2

plt.figure(figsize=(4,4))
#dust_plot(filenum,range(15,19),zoomangle=6.,zoomr=2)
#plt.savefig('output_plots/dust.pdf',bbox_inches='tight')

dust_densplot(filenum,2,zoomangle=1.,zoomr=6,maxcounts=40)
#plt.savefig('output_plots/dust.pdf',bbox_inches='tight')

# dr_corot = 1.05*np.sqrt(M_ratio/0.05)
# plt.axhline(y=1-dr_corot, color = 'darkgreen', linestyle = '--')
# plt.axhline(y=1+dr_corot, color = 'darkgreen', linestyle = '--')

# %% Plot autocorrelation

nice_plots()

t_stat = 200
size_bin = 19

time, torque_normalized, counts, N_bins, bins = read_force()

avg_counts = np.mean(counts[time>t_stat,:],axis=0)

Deltat=time[-1]/len(time)

plt.figure(figsize=(8,3))
plt.plot(time[-1500:-1],torque_normalized[-1500:-1,size_bin],markerfacecolor='white')
plt.grid()
plt.xlabel('time $t \ [2\pi / \\Omega_p]$')
plt.ylabel(r'normalized torque $ \tilde{\xi} \ ( t )$')
plt.axhline(y=0, color='k', linewidth=0.5)
plt.ylim(-0.1,0.1)
plt.xlim(time[-1500],time[-1])
#plt.savefig('output_plots/xi_zoom.pdf',bbox_inches='tight')

acf,tau = autocorrelation(torque_normalized[time>t_stat,size_bin])

plt.figure(figsize=(6,4))
plt.plot(time[0:3000],acf[0:3000], markerfacecolor='white')
plt.axvline(tau*Deltat, color='darkgreen', linestyle='--')
plt.text(tau*Deltat*1.2, 0.5, r"$t_{\ast} = $ "+f'{tau*Deltat:.3f}', color='darkgreen', verticalalignment='center')
plt.axhline(y=0, color='k', linewidth=0.5)
plt.grid(which='both')
plt.xlim(time[1],time[3000])
plt.xscale('log')
plt.xlabel('time lag $ k\Delta t \ [2\pi / \Omega_p]$')
plt.ylabel('torque autocorrelation $\\eta\ ( k )$')
#plt.savefig('output_plots/autocorrelation_time.pdf',bbox_inches='tight')

# %% gaussian distribution and error convergence

size_bin = 5
t_stat = 200
nslices = 100

plt.figure(figsize=(6,4))
counts, bins, _ = plt.hist(torque_normalized[time>t_stat,size_bin], bins=100, density=True, alpha=0.6, color='limegreen', label = r'$\tilde{\xi} \ $ counts')

avg, sigma, err = stats(torque_normalized[time>t_stat,size_bin])

bin_centers = 0.5 * (bins[1:] + bins[:-1])
pdf = norm.pdf(bin_centers, avg, sigma)

# Plot the fitted Gaussian
plt.plot(bin_centers, pdf, 'darkgreen', linestyle='--', label=f'Gaussian fit')
# plt.text(avg + 4*sigma, 0.5*max(pdf), r'std $( \ \tilde{\xi} \ ) =$'+f'${sigma:.2e}$',horizontalalignment='right')
# plt.text(avg + 4*sigma, 0.6*max(pdf), r'$\langle \ \tilde{\xi} \ \rangle=$'+f'${avg:.2e}$',horizontalalignment='right')
# plt.text(avg + 4*sigma, 0.4*max(pdf), r'$\delta \xi \ =$'+f'${3*err:.2e}$',horizontalalignment='right')
plt.axvline(x=0, color='k', linewidth=0.5)
plt.fill_betweenx([0, max(pdf)], avg - err, avg + err, color= 'pink', label=r'$ \left\langle \tilde{\xi} \right\rangle \pm \delta \xi $')
plt.plot([avg, avg], [0, norm.pdf(avg, avg, sigma)], color= 'firebrick',linestyle='-',label=r'$\left\langle \tilde{\xi} \right\rangle$')
plt.legend()
plt.xlabel(r"particle-averaged normalized torque $\tilde{\xi}$")
plt.ylabel("distribution density")
plt.grid()
#plt.savefig('output_plots/gaussian_distribution.pdf',bbox_inches='tight')

print(f"Average: {avg:.3e} +/- {err:.3e}")
print(f"STD: {sigma:.3e}")

NDt = np.zeros(nslices)
torqueavg = np.zeros(nslices)
torquestd = np.zeros(nslices)
torqueerr = np.zeros(nslices)
tsamples = np.zeros(nslices)

for i in range(nslices):

    tend= (time[-1]-t_stat)*(i+1)/(nslices) + t_stat
    #tstart= (time[-1]-tstat)*(i)/(nslices) + tstat

    time_filter = (time >= t_stat) & (time < tend)

    torqueavg[i], torquestd[i], torqueerr[i] = stats(torque_normalized[time_filter,size_bin])
    NDt[i] = len(torque_normalized[time_filter,size_bin])
    tsamples [i] = tend

plt.figure(figsize=(8,3))
plt.fill_between(tsamples, torqueavg - torqueerr, torqueavg + torqueerr, color='pink', label=r'$\left\langle \tilde{\xi} \right\rangle \pm \delta \xi$')
plt.plot(tsamples,torqueavg, label=r'$\left\langle \tilde{\xi} \right\rangle$')
plt.axhline(y=0, color='k', linewidth=0.5)
plt.legend(loc='lower right')
plt.grid()
plt.xlabel('time $t \ [2\pi / \\Omega_p]$')
plt.xlim(tsamples[0],tsamples[-1])
#plt.ylim(-0.015,0.015)
plt.ylabel('time-average of\n'+r'normalized torque $\left\langle \tilde{\xi}  \right\rangle$')
#plt.savefig('output_plots/error_convergence.pdf',bbox_inches='tight')

# %% Dust torque plot

fig,ax=plt.subplots(1,1,figsize=(6,3))
plt.axhline(y=0, color='k', linewidth=0.5)

output_dir = r'./problem/out_C9/'

time, torque_normalized, counts, N_bins, bins = read_force()

plot_sym_torque(Color="steelblue",Label='simulation')

# output_dir = r'./problem/out_C0.30_nd/'
# plot_sym_torque(Color="goldenrod",Label='$r_s=H_d$')
#stokes_ticks(ax)
mass_string = f'{M_ratio[0]/3.333e-6:.2f}'
plt.text(0.5, 0.02, f'$M_p ='+mass_string+'$ M$_\oplus$', ha='center', va='bottom',transform=ax.transAxes,)

plt.xscale('log')
plt.axhline(y=0, color='k', linewidth=0.7)
plt.grid(which='both')
#plt.xlabel('particle size $s$ [cm]')
#plt.ylabel(r'Dust torque $\Gamma_d/\Gamma_\ast$ ($\epsilon = 0.01$)')
plt.xlim(bins[0]*0.9,bins[-1]*1.1)
ax.set_yscale('symlog', linthresh=2)
ax.set_yticks([-1,-0.5,0,0.5,1,1.5,2,5,10,20,50])
ax.set_yticklabels(['-1','-0.5','0','0.5','1','1.5','2','5','10','20','50'])
minortickpos=[2,3,4,6,7,8,9,30,40]
# ax.set_yticks(minortickpos, minor=True)
# ax.set_yticks([-1,0,1,2,3,4,5,6,7,8,9,10,20,30,40,50])
# ax.set_yticklabels(['-1','0','1','2',' ',' ','5',' ',' ',' ',' ','10','20',' ',' ','50'])
plot_BLP18()
torque_g = torque_gas(-0.5)
plt.axhline(y=-torque_g/torque_0(), color='firebrick', linestyle='-.',label=r'$\Gamma_d+\Gamma_g=0$')
ax.set_ylim(bottom=-1)
#plt.legend(loc='upper left')

plt.savefig('output_plots/torque_Hd/torque'+mass_string+'.pdf',bbox_inches='tight')
# %% save data
folders=['out_C0.30',
         'out_C0.45',
         'out_C0.60',
         'out_C1',
         'out_C1.5',
         'out_C2.1',
         'out_C3',
         'out_C4',
         'out_C6',
         'out_C9']

# folders=['out_A0.30',
#          'out_A1',
#          'out_A3',
#          'out_A9']

filename = 'torques_Hill_radius.p'

Mass_ratio=np.zeros(len(folders))

for i in range(len(folders)):

    output_dir = r'./problem/'+folders[i]
    get_sym_data(printing=False)
    Mass_ratio[i] = M_ratio
    time, torque_normalized, counts, N_bins, bins = read_force()
    tstat = 100

    dust_gas_ratio = 0.01
 
    if i==0:
        torque=np.zeros((len(folders),N_bins))
        torque_err=np.zeros((len(folders),N_bins))

    for j in range(N_bins):

        torque_factor=torque_rescale_factor(bins[j],dust_gas_ratio)
        xi_avg, xi_sigma, xi_err = stats(torque_normalized[time>tstat,j]*torque_factor)
        torque[i,j]= xi_avg/torque_0()
        torque_err[i,j] = xi_err/torque_0()

# pk.dump((Mass_ratio, bins, torque, torque_err), open(f'output_data/'+filename, 'wb'))    

# Save torque data as CSV table

# Create header with size bins
header = ['Mass_ratio'] + [f'{b:.4e}' for b in bins]

# Create data rows with mass ratios and torque values
data = []
for i in range(len(folders)):
    row = [f'{Mass_ratio[i]:.4e}']
    for j in range(N_bins):
        row.append(f'{torque[i,j]:.4e}±{torque_err[i,j]:.4e}')
    data.append(row)

# Create DataFrame and save
df = pd.DataFrame(data, columns=header)
df.to_csv('output_data/torques_table.csv', index=False)
print(df)

# %% Size spectrum torque

filename = 'torques_dust_scaleheight.p'

Mass_ratios, size_bins, torque, torque_err = pk.load(open(f'output_data/'+filename, 'rb'))

def dust_torque_interp(size, Mass_ratio):
    
    #return torque and error for a given size and mass ratio

    #interpolate in log-log space
    log_M_ratios = np.log10(Mass_ratios)
    log_sizes = np.log10(size_bins)
    interp_func = interp2d(log_sizes, log_M_ratios, torque, kind='linear')
    interp_err_func = interp2d(log_sizes, log_M_ratios, torque_err, kind='linear')

    if size < size_bins[0]:
       torque_interp = 0.
       err_interp = 0.
    elif size > size_bins[-1]:
       torque_interp = 0.
       err_interp = 0.
       print("Warning: extrapolating torque above max size simulated")
    else:
        torque_interp = interp_func(np.log10(size), np.log10(Mass_ratio))
        err_interp = interp_err_func(np.log10(size), np.log10(Mass_ratio))

    return torque_interp, err_interp

smax = 1 #maximum size of distribution [cm]

def torque_distrib(Beta, Smax):

    dust_to_gas_ratio = 0.01

    int_slices = 200

    int_bin_boundaries = np.linspace(0, smax, int_slices+1)

    torque_integrated = np.zeros(len(Mass_ratios))
    err_integrated = np.zeros(len(Mass_ratios))

    for i in range(len(Mass_ratios)):

        for k in range(int_slices):

            size = 0.5*(int_bin_boundaries[k]+int_bin_boundaries[k+1])
            torque_interp, err_interp = dust_torque_interp(size, Mass_ratios[i])
            weight = (Beta+4)/Smax * (size/Smax)**(Beta+3) * (dust_to_gas_ratio/0.01)
            torque_integrated[i] += torque_interp * weight * (int_bin_boundaries[k+1]-int_bin_boundaries[k])
            err_integrated[i] += err_interp * weight * (int_bin_boundaries[k+1]-int_bin_boundaries[k])

    return torque_integrated, err_integrated

torque_frag, err_frag = torque_distrib(-3.5, smax)
torque_drift, err_drift = torque_distrib(-2.5, smax)

plt.figure(figsize=(6,4))
    
Mass_planet = Mass_ratios / 3.33e-6

plt.plot(Mass_planet, torque_frag,'--',color='navy',label=r'$\beta=-3.5$',linewidth=2)
plt.fill_between(Mass_planet, torque_frag - err_frag, torque_frag, alpha=0.3,color='gray',linewidth=0)
plt.plot(Mass_planet, torque_drift,'-',color='seagreen',label=r'$\beta=-2.5$',linewidth=2)
plt.fill_between(Mass_planet, torque_drift, torque_drift + err_drift, alpha=0.3,color='gray',linewidth=0)
plt.fill_between(Mass_planet, torque_frag, torque_drift, alpha=0.5,color='teal',linewidth=0)

plt.xscale('log')
plt.axhline(y=0, color='k', linewidth=0.5)
plt.axhline(y=-torque_gas(-0.5)/torque_0(),linestyle='-.',color='firebrick',label=r'$\Gamma_d+\Gamma_g=0$')
plt.xlabel(r'Planet mass $M_p$ [M$_\oplus$]')
plt.ylabel(r'Normalized dust torque $\Gamma_d/\Gamma_\ast$')
plt.grid(which='both')
plt.legend()

plt.text(0.5, 1.37, f'$s_c='+str(smax)+'$ cm', ha='center', va='top',transform=ax.transAxes,size=16)

plt.savefig('output_plots/torque_size_dist/torque'+str(smax)+'leg.pdf',bbox_inches='tight')



# %%
