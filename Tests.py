# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 18:31:32 2021

@author: Gerry O'Brien
"""
import GeneralBaseClass as g
import RayClass as r
import numpy as np
import matplotlib.pyplot as plt

# Testing propagation before fixing output
ray = r.Ray(pos=[1, 1, 0], dirn=[0, 0, 1])
prop_ray = g.SphericalRefraction()
prop_ray.propagate_ray(ray)
print(ray.p())
print(ray.k())
print(ray.vertices())

plt.xlabel('z (mm)')
plt.ylabel('y (mm)')
plt.grid()
plt.plot(np.array(ray.vertices())[:, 2], np.array(ray.vertices())[:, 1])

# %%
# Creating 11 paraxial rays to test focal length
# against spherical focal length formula before fixing output
for i in np.arange(-0.1, 0.12, 0.02):
    ray = r.Ray(pos=[0.0, i, 0.0], dirn=[0.0, 0.0, 1.0])
    prop_ray = g.SphericalRefraction()
    prop_ray.propagate_ray(ray)
    plt.xlabel('z (mm)')
    plt.ylabel('y (mm)')
    plt.grid()
    plt.plot(np.array(ray.vertices())[:, 2], np.array(ray.vertices())[:, 1])
    # plt.xlim(195, 202)  # checking the focal length by zooming in

# %%
# Testing focus from a point source - thin lens formula test
# Pretty bad approximation!
for i in np.arange(-0.1, 0.12, 0.02):
    ray = r.Ray(pos=[0.0, 0.0, 0.0], dirn=[0.0, i, 1.0])
    prop_ray = g.SphericalRefraction()
    prop_ray.propagate_ray(ray)
    output = g.OutputPlane(500)
    output.propagate_ray(ray)
    plt.title('Testing the thin lens formula')
    plt.xlabel('z (mm)')
    plt.ylabel('y (mm)')
    plt.grid()
    plt.plot(np.array(ray.vertices())[:, 2], np.array(ray.vertices())[:, 1])

# %%
# Refracting the bundle before fixing the ouput
for i in range(len(g.bun)):  # Plotting trace of bundle
    ray = r.Ray(pos=(g.bun[i]), dirn=[0.0, 0.0, 1.0])
    prop_ray = g.SphericalRefraction()
    prop_ray.propagate_ray(ray)
    # print(ray.vertices())
    plt.xlabel('z (mm)')
    plt.ylabel('y (mm)')
    plt.grid()
    plt.plot(np.array(ray.vertices())[:, 2], np.array(ray.vertices())[:, 1], color='purple')
    # plt.plot(ray.vertices()[2][0],ray.vertices()[2][1])
plt.show()

dot_positions = []
for i in range(len(g.bun)):  # Plotting bundle cross section at output
    ray = r.Ray(pos=(g.bun[i]), dirn=[0.0, 0.0, 1.0])
    prop_ray = g.SphericalRefraction()
    prop_ray.propagate_ray(ray)
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

# %%
# Testing propagation after fixing ouput plane


def parax_focus():  # Function finding the paraxial focal point
    parax_ray = r.Ray(pos=[0, 0.1, 0], dirn=[0, 0, 1])  # automating the focus
    refractor = g.SphericalRefraction()
    refractor.propagate_ray(parax_ray)
    l = -parax_ray.p()[1]/(parax_ray.k()[1]/np.linalg.norm(parax_ray.k()))
    z = parax_ray.p()[2]+l*(parax_ray.k()[2]/np.linalg.norm(parax_ray.k()))
    return z  # focal point


ray = r.Ray(pos=[1, 1, 0], dirn=[0, 0, 1])
refractor = g.SphericalRefraction()
refractor.propagate_ray(ray)
print(ray.p())
print(ray.k())
print(ray.vertices())
output = g.OutputPlane(parax_focus())
output.propagate_ray(ray)
print(ray.p())
print(ray.k())
print(ray.vertices())

plt.xlabel('z (mm)')
plt.ylabel('y (mm)')
plt.grid()
plt.plot(np.array(ray.vertices())[:, 2], np.array(ray.vertices())[:, 1])
# %%
# Testing no interception
ray = r.Ray(pos=[0, 50, 0], dirn=[0, 0, 1])
refractor = g.SphericalRefraction()
refractor.propagate_ray(ray)
# %%
# Refracting bundle after fixing the output plane

for i in range(len(g.bun)):  # Plotting trace of bundle
    ray = r.Ray(pos=(g.bun[i]), dirn=[0.0, 0.0, 1.0])
    refractor = g.SphericalRefraction()
    refractor.propagate_ray(ray)
    output = g.OutputPlane(parax_focus())
    output.propagate_ray(ray)
    # print(ray.vertices())
    plt.xlabel('z (mm)')
    plt.ylabel('y (mm)')
    plt.grid()
    plt.plot(np.array(ray.vertices())[:, 2], np.array(ray.vertices())[:, 1], color='purple')
    # plt.plot(ray.vertices()[2][0],ray.vertices()[2][1])
plt.show()

dot_positions = []
for i in range(len(g.bun)):  # Plotting bundle cross section at output
    ray = r.Ray(pos=(g.bun[i]), dirn=[0.0, 0.0, 1.0])
    refractor = g.SphericalRefraction()
    refractor.propagate_ray(ray)
    output = g.OutputPlane(parax_focus())
    output.propagate_ray(ray)
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
