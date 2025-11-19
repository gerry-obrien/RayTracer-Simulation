# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 15:45:15 2021

@author: Gerry O'Brien
"""

import GeneralBaseClass as g
import RayClass as r
import numpy as np
import matplotlib.pyplot as plt


# Generating a bundle and plotting the cross section at the input
bun = g.bundlefunc(2.5, 8)  # Radius set at 2.5mm
bun.append([0.0, 0.0, 0.0])
plt.title('Cross section of the ray bundle at the input')
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')
plt.plot(np.array(bun)[:, 0], np.array(bun)[:, 1], '.', color='purple')
plt.show()

def parax_focus():
    """
    Function to find the paraxial focal point of the spherical lens
    """
    parax_ray = r.Ray(pos=[0, 0.1, 0], dirn=[0, 0, 1])  # automating the focus
    refractor = g.SphericalRefraction()
    refractor.propagate_ray(parax_ray)
    l = -parax_ray.p()[1]/(parax_ray.k()[1]/np.linalg.norm(parax_ray.k()))
    z = parax_ray.p()[2]+l*(parax_ray.k()[2]/np.linalg.norm(parax_ray.k()))
    return z  # focal point


print('Paraxial focal point =', parax_focus())
# Plotting trace of bundle through the spherical lens
for i in range(len(bun)):
    ray = r.Ray(pos=(bun[i]), dirn=[0.0, 0.0, 1.0])
    refractor = g.SphericalRefraction()
    refractor.propagate_ray(ray)
    # Use the paraxial focal point as the ouput plane
    output = g.OutputPlane(parax_focus())
    output.propagate_ray(ray)
    # print(ray.vertices())
    plt.title('Ray trace for the ray bundle through the spherical refractor')
    plt.xlabel('z (mm)')
    plt.ylabel('y (mm)')
    plt.grid()
    plt.plot(np.array(ray.vertices())[:, 2], np.array(ray.vertices())[:, 1], color='purple')
    # plt.plot(ray.vertices()[2][0],ray.vertices()[2][1])
plt.show()

# Make a list to store the x, y coords of the rays at the ouput
dot_positions = []
# Plotting bundle cross section at output
for i in range(len(bun)):
    ray = r.Ray(pos=(bun[i]), dirn=[0.0, 0.0, 1.0])
    refractor = g.SphericalRefraction()
    refractor.propagate_ray(ray)
    output = g.OutputPlane(parax_focus())
    output.propagate_ray(ray)
    plt.title('Cross section of the ray bundle at the output')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    # print(ray.vertices())
    plt.plot(ray.p()[0], ray.p()[1], '.', color='purple')
    dot_positions.append([ray.p()[0], ray.p()[1]])
plt.show()

dot_positions = np.array(dot_positions)
dot_pos_squared = []
for vector in dot_positions:  # Finding RMS of output points
    r_squared = sum(vector**2)
    dot_pos_squared.append(r_squared)
rms = np.sqrt(np.mean(dot_pos_squared))
print('Root Mean Squared of points at output = ', rms, ' mm')
print('Diffraction scale of a laser with wavelength 588 nm = 0.0118 mm')
