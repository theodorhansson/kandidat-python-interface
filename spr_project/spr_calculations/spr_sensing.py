### Packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import os

from spr_project.spr_calculations.sensor_setup import plot_sensor_setup

### Defiinitions
PI = np.pi
DEG_TO_RAD = PI/180
RAD_TO_DEG = 180/PI
theta = np.arange(0,90,0.1)*DEG_TO_RAD
MM = 1e-3
UM = 1e-6
NM = 1e-9

### Plot arguments
fontsize_title  = 18
fontsize_label  = 16
fontsize_ticks  = 14
figure_width    = 10
figure_height   = 8
fontsize_legend = 12

lam0   = 984*NM
eps_Au = -40.650 + 1j*2.2254

def theta_spr(eps_metal, n_glass, n_analyte):
    sin_theta_spr = (1/n_glass)*np.sqrt((np.abs(eps_metal)*n_analyte**2)/(np.abs(eps_metal) - n_analyte**2))
    theta_spr = np.arcsin(sin_theta_spr)
    return theta_spr

def n_spr(eps_Au, n_prism, theta_spr):
    param = np.sin(theta_spr*np.pi/180)
    eps_analyte = (np.abs(eps_Au)*n_prism**2*param**2 )/(np.abs(eps_Au) + n_prism**2*param**2)
    return np.sqrt(eps_analyte)

def x_spr_detector(theta_spr, glass_thickness):
    return 2*glass_thickness*np.tan(theta_spr)

def spr_ang(x, glass_thickness):
    return np.arctan(x/(2*glass_thickness))

def n_glycerol(concentration, n_water):
    return concentration*1.474 + (1-concentration)*n_water

def pen_depth_intensity(lam, n_sub, n_bulk, theta_spr):
    return lam/(4*np.pi*np.sqrt(n_sub**2*np.sin(theta_spr)**2 - n_bulk**2))
    
def pen_depth(lam, n_sub, n_bulk, theta_spr):
    return lam/(2*np.pi*np.sqrt(n_sub**2*np.sin(theta_spr)**2 - n_bulk**2))

def angle_grating_equation(lam, n, period):
    return np.arcsin(lam/(n*period))


#%% Sensor properties

n_glass = 1.51
grating_period = 730*NM
glass_thickness = 1.4*MM

deflection_angle = angle_grating_equation(lam0, n_glass, grating_period)
angle_FWHM       = 15
deflection_from  = (deflection_angle*RAD_TO_DEG - angle_FWHM)*DEG_TO_RAD
deflection_to    = (deflection_angle*RAD_TO_DEG + angle_FWHM)*DEG_TO_RAD
print('--------------------------------')
print('Current design meant to center deflection at: ' + str(round(deflection_angle*RAD_TO_DEG, 1)))

n_low = 1.33
water_theta_spr_rad = theta_spr(eps_Au, n_glass, n_low)
water_theta_spr_deg = water_theta_spr_rad*RAD_TO_DEG
print('Predicted SPR for water: ' + str(round(water_theta_spr_deg, 2)))

n_high = 1.44
glycerol_theta_spr_rad = theta_spr(eps_Au, n_glass, n_high)

## Location of dip from vcsel row

water_dip_loc = x_spr_detector(water_theta_spr_rad, glass_thickness)
glycerol_dip_loc = x_spr_detector(glycerol_theta_spr_rad, glass_thickness)
print('--------------------------------')
print(f'Reference water dip will be offset -- ' + str(round(water_dip_loc/MM, 2)) + ' mm -- from vcsel row')
print('Highest possible detection dip will be offset -- ' + str(round(glycerol_dip_loc/MM, 2)) + ' mm -- from vcsel row')
print('--------------------------------')
print('Extent of detection: ' + str(round((np.abs(water_dip_loc - glycerol_dip_loc))/MM, 2)) + ' mm')
print('Extent of detection: ' + str(round(n_high - n_low, 2)) + ' RIU in bulk')


## Plot sensor setup
plot_sensor_setup(glass_thickness,
                  water_theta_spr_rad, glycerol_theta_spr_rad,
                  n_low, n_high,
                  deflection_angle,
                  deflection_to, deflection_from)

#%% Bulk sensitivity
## Sweep over refractive indicies
n_low  = 1.33
n_high = 1.40
n_sweep = np.linspace(n_low, n_high, 100)

## Corresponding angles
angle_sweep = theta_spr(eps_Au, n_glass, n_sweep)

## Different thickness of glass slides
glass_thickness = np.array([0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])*MM

loc_detector = np.zeros(shape=(len(glass_thickness), len(n_sweep)))
sensitivity_thickness = np.zeros(shape=(len(glass_thickness), len(n_sweep)))

fig, ax = plt.subplots(1,2, figsize=(figure_width, figure_height))
for i, thickness in enumerate(glass_thickness):
    
    loc_detector[i, :] = x_spr_detector(angle_sweep, thickness)
    sensitivity_thickness[i, :] = np.gradient(loc_detector[i, :], n_sweep)
    
    ax[0].plot(n_sweep, loc_detector[i, :]/MM,          label=f'Thickness {round(thickness/MM, 4)}')
    ax[1].plot(n_sweep, sensitivity_thickness[i, :]/MM, label=f'Thickness {round(thickness/MM, 4)}')

ax[0].set_xlabel('Refractive index [-]', fontsize=fontsize_label)
ax[1].set_xlabel('Refractive index [-]', fontsize=fontsize_label)
ax[0].set_ylabel('Postition on detector [mm]', fontsize=fontsize_label)
ax[1].set_ylabel('Sensitivity [mm/RIU]', fontsize=fontsize_label)

ax[0].legend(fontsize=fontsize_legend)
ax[0].set_title('Detector Position', fontsize=fontsize_title)
ax[1].legend(fontsize=fontsize_legend)
ax[1].set_title('Sensitivity', fontsize=fontsize_title)

plt.tight_layout()

#%% Penetration depth in bulk

n_analyte_low  = 1.33
n_analyte_high = 1.44
n_analyte_sweep = np.linspace(n_analyte_low, n_analyte_high, 100)

## Corresponding angles
angle_sweep = theta_spr(eps_Au, n_glass, n_analyte_sweep)
pen_depth_sweep = pen_depth(lam0, n_glass, n_analyte_sweep, angle_sweep)

fig, ax = plt.subplots(1,1)
ax.plot(n_analyte_sweep, pen_depth_sweep/NM)

ax.set_xlabel('Refractive index [-]', fontsize=fontsize_label)
ax.set_ylabel('Penetration depth [nm]', fontsize=fontsize_label)

ax.legend(fontsize=fontsize_legend)
ax.set_title('Penetration depth', fontsize=fontsize_title)
plt.tight_layout()







