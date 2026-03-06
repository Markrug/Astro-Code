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

def makestar(center_density, gridsize):
    
    n_e = center_density/(2*m_h)
    P_gas = K * n_e**(5/3)
    mass = 0
    radius = gridsize

    while P_gas > 0:
        
        n_e = (P_gas/K)**(3/5)
        density = 2*m_h*n_e
        
        #conservation mass
        dmdr = 4*pi*radius**2 * density
        #hydrostatic equilibrium
        dPdr = -density * ((G * mass) / (radius**2))
        
        mass = mass + dmdr*gridsize
        P_gas = P_gas + dPdr*gridsize
        radius = radius + gridsize
        
    return mass, radius


m_sun = M_sun.to(u.kg).value
r_sun = R_sun.to(u.m).value
Masses = []
Radii = []

#loop through central densities until 0.2 solar mass found, save 0.2 - 1.4 values
central_densities = np.logspace(7, 17, num=1000)
for density in central_densities:
    mass, radius = makestar(density, 500)

    mass_solar = mass/m_sun

    if mass_solar >= 0.2:
        Masses.append(mass_solar)
        Radii.append(radius/r_sun)

    if mass_solar > 1.4:
        break

#convert radii to 1e3 km
Radii_1000km = np.array(Radii)*r_sun/1e6

#plot mass-radius relation
plt.plot(Radii_1000km, Masses)
plt.ylabel('Mass (Solar Masses)')
plt.xlabel('Radius (1000 km)')
plt.title('Mass-Radius Relation for White Dwarfs')
plt.grid()
plt.show()







