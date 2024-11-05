###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
#                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
#                    Jasper Leonora Kamermans (j.l.p.d.kamermans@students.uu.nl)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################

# Import packages

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
# matplotlib.use("Agg")
from pylab import *
import h5py
import numpy as np
import glob
import random
# from iminuit import Minuit
# from iminuit.cost import ExtendedUnbinnedNLL, BinnedNLL, UnbinnedNLL
import scipy.odr as odr
# import scipy.integrate as integrate
# from jacobi import jacobi
np.random.seed(1)

np.seterr(divide='ignore', invalid='ignore') #  THIS MAKES IT SO THE ERRORS DO NOT SHOW

# The function of the spike
def logFit(B,x):
    return (x/a/B[0])**B[3] + B[2]

def fhernq(a, x_data):
    return Mhalo/(2*np.pi)*a/x_data/(x_data+a)**3

# Plot parameters
params = {
    "font.size": 20,
    "figure.figsize": (8, 8),
    "figure.subplot.left": 0.0,
    "figure.subplot.right": 1,
    "figure.subplot.bottom": 0,
    "figure.subplot.top": 1,
    "figure.subplot.wspace": 0,
    "figure.subplot.hspace": 0,
    "lines.markersize": 6,
    "lines.linewidth": 3.0,
    "figure.max_open_warning": 0,
    "mathtext.fontset": 'cm' 
}
rcParams.update(params)

# Make it so we can use the correct fonts
from matplotlib import font_manager

font_dir = ['/usr/share/texlive/texmf-dist/fonts/opentype/']
for font in font_manager.findSystemFonts(font_dir):
    try:
        font_manager.fontManager.addfont(font)
    except:
        continue
   
plt.rcParams["font.family"] = "TeX Gyre Termes"
plt.clf()

# The binning of the data!
binsize = .05
begin_bins = -6
end_bins = 1

# Radial bins
radial_bins = np.arange(begin_bins, end_bins, binsize)

# DETERMINE THE NUMBER OF SNAPSHOTS
newest_snap_name = max(glob.glob('box_*.hdf5'))#, key=os.path.getctime)
n_snapshots = int(newest_snap_name.replace('box_','').replace('.hdf5','')) + 1


# DETERMINE THE START TIMESTEP FOR EACH SNAPSHOT
dt = .01
# Spike radius is stored here
SpikeRadius = []

# Uncertainty in the spike slope
SigmaSlope = []

# Save the fitted parameters

list1 = []
list2 = []
list3 = []
list4 = []
mulist = []

tempie = []
# empty array for the fitted parameters
empty = np.array([0,0,0,0])

# Go over every datafile
for j in range(n_snapshots):
        print(j)
        # Try and read the datafile
        try:
            sim = h5py.File("box_%04d.hdf5"%j, "r")
            # Read in the softening (should be constant over all particles)
            softening = sim["/PartType1/Softenings"][:][0]
            # Read the data of the particles from the file
            DMpos = sim["/PartType1/Coordinates"][:]
            DMvel = sim["/PartType1/Velocities"][:,:]
            massParticle = sim["/PartType1/Masses"][:][0]*1e10
            BHpos = sim["/PartType5/Coordinates"][:]
            Mhalo = len(DMpos)*massParticle
            
            # Calculate the relative location of the DM particels with respect to the Black Hole
            # The commented part calculates it relative to the Centre Of Mass of the Halo itself
            DMposrel = DMpos - BHpos
            BHvel = sim["/PartType5/Velocities"][:]
            massBH = sim["/PartType5/DynamicalMasses"][:][0]*1e10
         
        # If it fails, it is probably corrupted!
        except:
            print("File "+ str(j) + " is corrupted!")
            continue
        
        
        # The radius from the BH for every DM particle
        r = np.sqrt((DMpos[:,0] - BHpos[:,0])**2 + (DMpos[:,1] - BHpos[:,1])**2 + (DMpos[:,2] - BHpos[:,2])**2)
           
        
        # Calculate mu!
        mu = massBH/(Mhalo + massBH) 
        BhHalo =  massBH/Mhalo
        print("mu = "+ str(mu))
        
        # Bin the radial data!
        bins, _ = np.histogram(r, bins = 10**radial_bins)       
        
        x_data = 10**radial_bins[:-1] + (10**radial_bins[1:] - 10**radial_bins[:-1])/2
        
        # Calculate the radial density if a single particle is present in each bin
        singleParticleDensity = massParticle/((4/3*np.pi*(10**radial_bins[1:])**3) - (4/3*np.pi*(10**radial_bins[:-1])**3))

        
        # Scale the data with the mass of the particle and the volume in each bin
        y_data = bins*singleParticleDensity
        y_data_std = np.sqrt(bins)*singleParticleDensity
        
        
        
        
        
        if j == 0:
            # Determine a

            rho_crit = 2.775e11 * 0.6777**2 #Msun/Mpc^3
            Rvir3 = Mhalo / ( (4. * np.pi / 3.) * 200 * rho_crit )
            Rvir = Rvir3 ** (1./3.) # Virial radius [Mpc units]
            Rvir *= 1e3 # [kpc units]

            z = 0
            alpha = 1.7543 - 0.2766*(1 + z) + 0.02039*(1+z**2)
            beta = 0.2753 + 0.00351*(1 + z) - 0.3038*(1 + z)**0.0269
            gamma = -0.01537 + 0.02102*(1 + z)**(-0.1475)

            c = 10**(alpha + beta*np.log10(Mhalo)*(1+ gamma*np.log10(Mhalo)**2))

            a = Rvir/c #kpc - scale radius
            
            continue
            
        # After the first datafile, we fit the spike!
        else:
            # Select the data that is below 25 times a
            y_datacool = y_data
            y_data_stdcool = y_data_std
            binscool = bins
            x_datacool = x_data
            
            # CHANGE THESE VALUES
            # The lower fitting range is declared here, you have to manualy mask 
            # away the shockwaves I'm afraid
            if j < 15:
                fittingradius = 0.1
            elif j < 25:
                fittingradius = 0.2
            elif j < 40:
                fittingradius = 0.3
            elif j < 40:
                fittingradius = 0.2
            
            lowerfittingradius = 2.5*softening
                
            # Mask away the data outside of the fitting range
            y_data = y_data[x_data < fittingradius]
            y_data_std = y_data_std[x_data < fittingradius]
            bins = bins[x_data < fittingradius]
            x_data = x_data[x_data < fittingradius]

            y_data = y_data[x_data > lowerfittingradius]
            y_data_std = y_data_std[x_data > lowerfittingradius]
            bins = bins[x_data > lowerfittingradius]
            x_data = x_data[x_data > lowerfittingradius]
            
            print("fitting bins = " + str(len(x_data)))

            # Calculate spike divided by Hernquist
            SpikeOverHernquist = y_data/fhernq(a, x_data)
            SpikeOverHernquist_std = y_data_std/fhernq(a, x_data)
            
            # Select model
            odr_model = odr.Model(logFit)
            
            # Fit a Spike profile on the data
            odr_data = odr.RealData(x_data, SpikeOverHernquist, sy = SpikeOverHernquist_std)
            # Initial values of the fit
            par_best_1 = [0,1,1,-1.333]

            # Perform fit!
            odr_obj = odr.ODR(odr_data, odr_model, beta0 = par_best_1, maxit = 10000)
            odr_obj.set_job(fit_type = 2)
            odr_res = odr_obj.run()
            
            # The results of the fit
            par_best = odr_res.beta
            par_sig_ext = odr_res.sd_beta
            par_cov_beta = odr_res.cov_beta
            par_best_1 = par_best.copy()
            chi = odr_res.res_var
            
            #If a spike was confidently found (change the number for more/less confidence)
            if logFit(par_best,x_data[0]) > 1.15:
                if  True:
                    # Add the sigma of every parameter to empty
                    empty = np.vstack((empty,par_sig_ext))
                    
                    # Append the best fitting values of the parameters
                    list1.append(par_best_1[0])
                    list2.append(par_best_1[1])
                    list3.append(par_best_1[2])
                    list4.append(par_best_1[3])
                    mulist.append(mu)
                        
                
            # Plot the data and fit
            plt.plot(x_datacool,y_datacool/fhernq(a, x_datacool), label = "Data", c = "navy")
            plt.fill_between(x_datacool, y_datacool/fhernq(a, x_datacool) - 2*y_data_stdcool/fhernq(a, x_datacool), y_datacool/fhernq(a, x_datacool) + 2*y_data_stdcool/fhernq(a, x_datacool), alpha = 0.3, color = "navy", label = r"$\mathdefault{2\sigma}$ uncertainty")
            plt.plot(x_datacool,logFit(par_best_1,x_datacool), label = "Best fit", c = "r", linestyle = "dashed")
            plt.yscale('log')
            plt.xscale('log')
            plt.vlines(lowerfittingradius, 1e-10, 1e20, 'purple', label = "Fitting radii", linestyle = "dashdot")
            plt.vlines(fittingradius, 1e-10, 1e20, 'purple',linestyle = "dashdot")
            plt.legend()
            plt.xlim(0.0005,1)
            plt.ylim(.8,10)
            plt.grid(True, which = "major", alpha = 0.4)
            plt.grid(True, which = "minor", ls = '--', alpha = 0.2)
            plt.xlabel('Radius [kpc]')
            plt.ylabel(r'$\mathdefault{\rho/\rho_{hernq}}$')
            plt.yticks([0.8,0.9,1,2,3,4,5,6,7,8,9,10], labels = ['0.8','','1','2','','4','','6','','8','','10'])
            plt.tight_layout(pad = 0.31)
            plt.title(j)
            plt.savefig("figure/fit" + str(j) + ".pdf")
            plt.show()
 

# Plot the best fits of the parameters
plt.plot(mulist[:],list1[:])
plt.title("B1") 
plt.show()
plt.plot(mulist[:],list3[:])
plt.title("B3")
plt.show()
plt.plot(mulist[:],1 - np.array(list4[:]))
# Plot the error values of the slope as well
plt.fill_between(mulist[:], 1 - np.array(list4[:]) - 2*np.array(SigmaSlope[:]), 1 - np.array(list4[:]) + 2*np.array(SigmaSlope[:]), alpha = 0.5)
plt.hlines(7/3,0,.1,'r', label = r"Gondola et Silk spike ($\gamma = \beta - 1 = -7/3$)")
plt.title("Value for slope parameter Beta")
plt.xlabel(r"$\mu$")
plt.ylabel(r"$\beta$")
plt.show()



