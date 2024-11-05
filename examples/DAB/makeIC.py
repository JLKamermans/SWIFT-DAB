###############################################################################
 # This file is part of SWIFT.
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

import h5py
import numpy as np
from scipy import integrate


np.random.seed(1)

binsize = .1
begin_bins = -5
end_bins = 0

radial_bins = np.arange(begin_bins, end_bins, binsize)
    
pi = np.pi

def func(E, haloMass, a, G):
    # x specific energy
    G = 4.30091e-6 # km^2/s^2 kpc Msun^-1
    q = np.sqrt((a * E)/(G * haloMass))
    vg = np.sqrt(G * haloMass / a)
    f = haloMass / (8. * np.sqrt(2.) * np.pi**3 * a**3 * vg**3 * (1.-q**2)**(5./2.))
    f *= (3. * np.arcsin(q) + q * np.sqrt(1.-q**2) * (1. - 2*q**2) * (8. * q**4 - 8. * q**2 - 3.))
    return f

# Generates a swift IC file #

G = 4.30091e-6 # km^2/s^2 kpc Msun^-1
boxsize = 1000 # some size
center = boxsize / 2.
numPart = int(100**3) # The number of particles in the simulation! We put it slightly lower for convenience
L = numPart**(1./3.)
eta = 1.2348          # 48 ngbs with cubic spline kernel


# Here the initial masses are given for the halo and the central BH, and the final BH mass
haloMass = 1e4 #Msun
part_mass = haloMass/numPart #Msun
bhMass = part_mass #Msun. If the initial BH mass is too large, the function from appendix A1 should be implemented instead.
mass_ratio = bhMass/(haloMass + bhMass)
finalMassBH = 1e3 #Msun
bhMassGrowthGyr = 1e3/4

# Now, calculate the scale mass based on the results of Correa et al. https://arxiv.org/abs/1502.00391
rho_crit = 2.775e11 * 0.6777**2 #Msun/Mpc^3
Rvir3 = haloMass / ( (4. * np.pi / 3.) * 200 * rho_crit )
Rvir = Rvir3 ** (1./3.) # Virial radius [Mpc units]
Rvir *= 1e3 # [kpc units]

z = 0
alpha = 1.7543 - 0.2766*(1 + z) + 0.02039*(1+z**2)
beta = 0.2753 + 0.00351*(1 + z) - 0.3038*(1 + z)**0.0269
gamma = -0.01537 + 0.02102*(1 + z)**(-0.1475)

c = 10**(alpha + beta*np.log10(haloMass)*(1+ gamma*np.log10(haloMass)**2))

a = Rvir/c #kpc - scale radius
print(c, a, Rvir)

# Set units
unit_length_in_cgs = 3.085678e21
unit_mass_in_cgs = 1.98848e43
unit_time_in_cgs = 3.085678e16
unit_current_in_cgs = 1
unit_temperature_in_cgs = 1

# Generate the radial coordinates of the particles
x = np.random.uniform(0, 0.98, numPart) 
r = a / ( np.sqrt(1./x)-1. )

# Output File
fileName = "test.hdf5"
file = h5py.File(fileName, 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = boxsize
grp.attrs["NumPart_Total"] =  [0, numPart, 0, 0, 0, 1]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [0, numPart, 0, 0, 0, 1]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0, part_mass / 1e10, 0.0, 0.0, 0.0, bhMass/1e10]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 3

#Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = unit_length_in_cgs
grp.attrs["Unit mass in cgs (U_M)"] = unit_mass_in_cgs
grp.attrs["Unit time in cgs (U_t)"] = unit_time_in_cgs
grp.attrs["Unit current in cgs (U_I)"] = unit_current_in_cgs
grp.attrs["Unit temperature in cgs (U_T)"] = unit_temperature_in_cgs

# Create data

# masses
mass = np.ones(numPart) * part_mass

# positions
theta = np.random.uniform(0, 1, numPart) * 2. * np.pi # uniform [0,2pi)
u = 2. * np.random.uniform(0, 1, numPart) - 1 #uniform [-1,1)


pos = np.zeros((numPart,3))
pos[:,0] =  r * np.cos(theta) * np.sqrt(1. - u**2)
pos[:,1] =  r * np.sin(theta) * np.sqrt(1. - u**2)
pos[:,2] =  r * u
# centering
pos[:,0] += center
pos[:,1] += center
pos[:,2] += center
radius = np.sqrt((pos[:,0]-center)**2 + (pos[:,1]-center)**2 + (pos[:,2]-center)**2)



P = G * haloMass /(r + a) + G*bhMass/r #Potential energy

# Generate distribution function
E = np.arange( -5, np.log10(np.max(P)), 0.01)
E = 10**E
f = func(E,haloMass,a,G)

# Make sure no nans
E = E[~np.isnan(f)]
f = f[~np.isnan(f)]



# Choose E from distribution function
n = numPart
Y = np.random.uniform(0, 1, numPart)
ans = np.ones(n)
for i in range(numPart):
    newE = E[E < P[i]] #newE goes to maxE(R)=P(R)
    F = integrate.cumtrapz(f[E<P[i]] * np.sqrt(P[i]-newE),newE, initial=0) #cummulative distribution
    newE = newE[~np.isnan(F)]
    F = F[~np.isnan(F)]
    ans[i] = np.interp(Y[i]*F[-1], F, newE) #find at which energy F=rand


abs_random_v = np.sqrt( 2. * (P - ans)) #velocity magnitude

theta = np.random.uniform(0, 1, numPart) * 2. * np.pi # uniform [0,2pi)
u = 2. * np.random.uniform(0, 1, numPart) - 1 #uniform [-1,1)

vel = np.zeros((numPart,3))
vel[:,0] =  abs_random_v * np.cos(theta) * np.sqrt(1. - u**2)
vel[:,1] =  abs_random_v * np.sin(theta) * np.sqrt(1. - u**2)
vel[:,2] =  abs_random_v * u

# Tset the BH to be at the centre and with no initial velocity
bhPos = [center,center,center]
bhVel = [0,0,0]


# Particle group
grp5 = file.create_group("/PartType5")
ds5 = grp5.create_dataset('Coordinates', (1, 3), 'd', data=bhPos)
ds5 = grp5.create_dataset('Velocities', (1, 3), 'f', data=bhVel)
ds5 = grp5.create_dataset('Masses', (1,1), 'f', data=bhMass/1e10)
ds5 = grp5.create_dataset('MassGrowthGyr', (1,1), 'f', data=bhMassGrowthGyr)
ds5 = grp5.create_dataset('MaxMass', (1,1), 'f', data=finalMassBH)


l = (4. * np.pi * radius**2 / numPart)**(1./3.) #local mean inter-particle separation
ds5 = grp5.create_dataset('SmoothingLength', (1,1), 'f', data=np.mean(l)/25)

ids = np.linspace(0, 1, 1, endpoint=False).reshape((1,1))
ds5 = grp5.create_dataset('ParticleIDs', (1, 1), 'L')
ds5[()] = ids + 1


# Particle group
grp = file.create_group("/PartType1")
ds = grp.create_dataset('Coordinates', (numPart, 3), 'd', data=pos)
ds = grp.create_dataset('Velocities', (numPart, 3), 'f', data=vel)
ds = grp.create_dataset('Masses', (numPart,1), 'f', data=mass/1e10)

l = (4. * np.pi * radius**2 / numPart)**(1./3.) #local mean inter-particle separation
print("one-twenty-fifth of mean inter partcle seperation: " + str(np.average(l)/25))
h = np.full((numPart, ), 1/25*l)
ds = grp.create_dataset('SmoothingLength', (numPart,1), 'f', data=h)

ids = np.linspace(0, numPart, numPart, endpoint=False).reshape((numPart,1))
ds = grp.create_dataset('ParticleIDs', (numPart, 1), 'L')
ds[()] = ids + 1

file.close()

N200 = len(np.ma.compressed(r[r < Rvir]))

# Calculate the best softening based on Power et al. The prefactor of 6 is slightly "try and find out"-based I'm afraid
print("The best softening length would be: " + str(6*Rvir/np.sqrt(N200)))
