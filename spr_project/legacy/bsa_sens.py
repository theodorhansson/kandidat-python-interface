#%% --- CALCULATION FOR BSA MEASUREMENT
import numpy as np
import matplotlib.pyplot as plt

DEG = np.pi/180

def ResonantAngle(eps_3r, n_glass, n_analyte):
    param = np.sqrt((n_analyte**2*np.abs(eps_3r))/(np.abs(eps_3r) - n_analyte**2))/n_glass
    theta_spr = np.arcsin(param)
    return theta_spr

def pen_depth(lam0, theta, n_sub, n_channel):
    return lam0/(2*np.pi*np.sqrt(n_sub**2*np.sin(theta)**2 - n_channel**2))

# Vacuum wavelength 
lam0 = 984

# Bulk sensitivity
S = 18945 #um/RIU

# Refractive index of proteins
n_protein = 1.45

# Refractive index of water
n_water = 1.33

# Refractive index of glass, substrate or prism
n_glass = 1.51

# Permativitty of gold
eps_Au = -40.650+ 1j*2.2254

n_Au = np.sqrt((1/2)*(np.sqrt(np.real(eps_Au)**2 + np.imag(eps_Au)**2) + np.real(eps_Au)))
k_Au = np.sqrt((1/2)*(np.sqrt(np.real(eps_Au)**2 + np.imag(eps_Au)**2) - np.real(eps_Au)))


# Thickness of a single protein layer
d = 7 #nm

# How many layers to calculate for
number_of_layers = np.arange(50)+1

# Initialize shift and effective refractive index arrays
dx        = np.zeros_like(number_of_layers, dtype=float)
neff      = np.zeros_like(number_of_layers, dtype=float) 
theta_spr = np.zeros_like(number_of_layers, dtype=float) 
ldep      = np.zeros_like(number_of_layers, dtype=float)

# SPR Angle for water
theta_spr[0] = ResonantAngle(eps_Au, n_glass, n_water)

# Penetration depth
ldep[0] = pen_depth(lam0, theta_spr[0], n_glass, n_water)


#%%
### STANDARD FORUMULA
# for i in range(len(number_of_layers)):
for i in range(0,3):
    # Calculate the SPR shift induced by the layer
    dx[i] = S*(n_protein - n_water)*(1 - np.exp(-2*d*(number_of_layers[i] - 1)/ldep[i]))

    # Calculate the **effective** refractive index within the field
    neff[i] = n_water + (n_protein - n_water)*(1 - np.exp(-2*d*(number_of_layers[i] - 1)/ldep[i]))
    
    # Calculate the new SPR angle with this new refractive index
    theta_spr[i] = ResonantAngle(eps_Au, n_glass, neff[i])

    # Calculate the new penetration depth with the new angle
    ldep[i] = pen_depth(lam0, theta_spr[i], n_glass, neff[i])

# Plot shift for each layer

# plt.plot(number_of_layers, theta_spr)1
plt.plot(number_of_layers, dx)

# plt.plot(m[0:-1], np.diff(dx))
# ax8.scatter(m,dx, edgecolor='b',facecolor='none') 
# ax8.set_xlabel('Layer No')
# ax8.set_ylabel('SPR Coordinate on Detector, um')