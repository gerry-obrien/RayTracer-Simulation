# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 14:34:56 2021

@author: Gerry O'Brien
"""

import GeneralBaseClass as g
import RayClass as r
import numpy as np
import matplotlib.pyplot as plt

# Generating a bundle and plotting the cross section at the input
bun = g.bundlefunc(5, 8)
bun.append([0.0, 0.0, 0.0])
plt.title('Cross section of the ray bundle at the input')
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')
plt.plot(np.array(bun)[:, 0], np.array(bun)[:, 1], '.', color='purple')
plt.show()


def parax_focus2():
    """
    Function to find the paraxial focal length
    of the plano-convex lens with the curved side towards the input
    """
    parax_ray = r.Ray(pos=[0, 0.1, 0], dirn=[0, 0, 1])
    convex = g.SphericalRefraction(z0=100, curvature=0.02, n1=1.0, n2=1.5168, aperture_radius=10)
    convex.propagate_ray(parax_ray)
    plano = g.SphericalRefraction(z0=105, curvature=0.0, n1=1.5168, n2=1.0, aperture_radius=10)
    plano.propagate_ray(parax_ray)
    l = -parax_ray.p()[1]/(parax_ray.k()[1]/np.linalg.norm(parax_ray.k()))
    z = parax_ray.p()[2]+l*(parax_ray.k()[2]/np.linalg.norm(parax_ray.k()))
    return z  # focal point


print('Paraxial focal point =', parax_focus2())
# Plotting the trace of the bundle through the plano-convex lens
# with the curved side towards the input
for i in range(len(bun)):
    ray = r.Ray(pos=(bun[i]), dirn=[0.0, 0.0, 1.0])
    convex = g.SphericalRefraction(z0=100, curvature=0.02, n1=1.0, n2=1.5168, aperture_radius=10)
    convex.propagate_ray(ray)
    plano = g.SphericalRefraction(z0=105, curvature=0.0, n1=1.5168, n2=1.0, aperture_radius=10)
    plano.propagate_ray(ray)
    output = g.OutputPlane(parax_focus2())
#    output = OutputPlane(400)
    output.propagate_ray(ray)
    # print(ray.vertices())
    plt.title('Ray trace through the Plano-Convex lens with the curved side towards the input')
    plt.xlabel('z (mm)')
    plt.ylabel('y (mm)')
    plt.grid()
    plt.plot(np.array(ray.vertices())[:, 2], np.array(ray.vertices())[:, 1], color='purple')
    # plt.plot(ray.vertices()[2][0],ray.vertices()[2][1])
plt.show()

dot_positions = []
# Plotting the bundle cross section at output for plano-convexlens
# with the flat side towards the input
for i in range(len(bun)):
    ray = r.Ray(pos=(bun[i]), dirn=[0.0, 0.0, 1.0])
    convex = g.SphericalRefraction(z0=100, curvature=0.02, n1=1.0, n2=1.5168, aperture_radius=10)
    convex.propagate_ray(ray)
    plano = g.SphericalRefraction(z0=105, curvature=0.0, n1=1.5168, n2=1.0, aperture_radius=10)
    plano.propagate_ray(ray)
#    output = g.OutputPlane(g.parax_focus())
    output = g.OutputPlane(parax_focus2())
    output.propagate_ray(ray)
    plt.title('Cross section at the output for the Plano-Convex lens with the curved side towards the input')
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
print('Diffraction scale of a laser with wavelength 588 nm = 0.00588 mm')
