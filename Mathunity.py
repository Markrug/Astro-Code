import math
import numpy as np
from scipy import constants
import astropy.units as u
from astropy.constants import M_sun, R_sun
import matplotlib.pyplot as plt


pi = math.pi
h = constants.h
G = constants.G
m_h = constants.m_p    #Assuming fully ionized hydrogen
m_e = constants.m_e

K = (h**2/(5*m_e)) * (3/(8*pi))**(2/3)

#scalars for order unity: values are 1e-1 - 1e1
M0 = M_sun.to(u.kg).value
R0 = 1e7
density0 = 1e8
n_e0 = density0/(2*m_h)
P_gas0 = K * n_e0**(5/3)

C1 = 4*pi*R0**3*density0/M0
C2 = G*M0*density0/(P_gas0*R0)

def makestar(center_density, gridsize):
    
    n_e = center_density * density0/(2*m_h)
    P_gas = (K * n_e**(5/3))/P_gas0
    mass = 0
    radius = gridsize

    while P_gas > 0:
        
        n_e = (P_gas*P_gas0/K)**(3/5)
        density = (2*m_h*n_e)/density0
        
        #conservation mass
        dmdr = C1 * radius**2 * density
        #hydrostatic equilibrium
        dPdr = C2 * (-density * mass) / (radius**2)
        
        mass = mass + dmdr*gridsize
        P_gas = P_gas + dPdr*gridsize
        radius = radius + gridsize
        
    return mass, radius


Masses = []
Radii = []
kept_central_densities = []

#loop through central densities until 0.2 solar mass found, save 0.2 - 1.4 values
central_densities = np.logspace(-1, 5, num=1000)
for density in central_densities:
    mass, radius = makestar(density, 0.001)

    if mass >= 0.2:
        Masses.append(mass)
        Radii.append(radius)
        kept_central_densities.append(density*density0/1000)

    if mass> 1.4:
        break

#convert radii to 1e3 km
Radii_1000km = np.array(Radii)*R0/1e6

#plot mass-radius relation and plot y = x^(1/3) for comparison as well
x=[]
for mass in Masses:
    x.append(mass**(-1/3))
    
scale = Radii_1000km[0]/x[0]
x = (np.array(x)+0.1)*scale

plt.plot(Radii_1000km, Masses)
plt.plot(x, Masses, linestyle="--", c='red')
plt.ylabel('Mass (Solar Masses)')
plt.xlabel('Radius (1000 km)')
plt.legend(["Computed Model", "M^(-1/3) Relation (offset added for clarity)"])
plt.title('Mass-Radius Relation for White Dwarfs')
plt.grid()
plt.show()


# plot core density vs mass
# core density quickly exceeds the relativistic limit, but that only holds for the center of the white dwarf
plt.plot(Masses, kept_central_densities)
plt.xlabel('Mass (Solar Masses)')
plt.ylabel('Central Density (g/cm^3)')
plt.title('Core Density vs Mass')
plt.axhline(1e6, linestyle="--", c='red')
plt.legend(["Computed Model", "Ultra-Relativistic Limit"])
plt.grid()
plt.show()

avg_densities=  [(Masses[i] * M0) / (4 * np.pi * (Radii[i] * R0)**3) * 1e-3 for i in range(len(Masses))]

# plot average density vs mass, stays well below limit as expected
plt.plot(Masses, avg_densities)
plt.xlabel('Mass (Solar Masses)')
plt.ylabel('Average Density (g/cm^3)')
plt.title('Average Density vs Mass')
plt.axhline(1e6, linestyle="--", c='red')
plt.legend(["Computed Model", "Ultra-Relativistic Limit"])
plt.grid()
plt.show()