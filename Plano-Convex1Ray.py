# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 13:39:23 2021

@author: Gerry O'Brien
"""

import GeneralBaseClass as g
import RayClass as r
import numpy as np
import matplotlib.pyplot as plt

# Tracing 1 ray through the plano-convex lens
ray = r.Ray(pos=[0, 2, 0], dirn=[0.0, 0.0, 1])

plano = g.SphericalRefraction(z0=100, curvature=0.0, n1=1.0, n2=1.5168, aperture_radius=10)
plano.propagate_ray(ray)
convex = g.SphericalRefraction(z0=105, curvature=-0.02, n1=1.5168, n2=1.0, aperture_radius=10)
convex.propagate_ray(ray)
print(ray.p(), ray.k())
output = g.OutputPlane(300)
output.propagate_ray(ray)

plt.xlabel('z (mm)')
plt.ylabel('y (mm)')
plt.grid()
plt.plot(np.array(ray.vertices())[:, 2], np.array(ray.vertices())[:, 1])

# %%
# Tracing 1 ray through the plano-convex lens reversed
ray = r.Ray(pos=[0, 2, 0], dirn=[0.0, 0.0, 1])

convex = g.SphericalRefraction(z0=100, curvature=0.02, n1=1.0, n2=1.5168, aperture_radius=10)
convex.propagate_ray(ray)
plano = g.SphericalRefraction(z0=105, curvature=0.0, n1=1.5168, n2=1.0, aperture_radius=10)
plano.propagate_ray(ray)
output = g.OutputPlane(300)
output.propagate_ray(ray)

plt.xlabel('z (mm)')
plt.ylabel('y (mm)')
plt.grid()
plt.plot(np.array(ray.vertices())[:, 2], np.array(ray.vertices())[:, 1])
# %%
# Tracing a few rays through the plano-convex lens
for i in np.arange(-0.1, 0.12, 0.02):
    ray = r.Ray(pos=[0.0, i, 0.0], dirn=[0.0, 0.0, 1.0])
    plano = g.SphericalRefraction(z0=100, curvature=0.0, n1=1.0, n2=1.5168, aperture_radius=10)
    plano.propagate_ray(ray)
    convex = g.SphericalRefraction(z0=105, curvature=-0.02, n1=1.5168, n2=1.0, aperture_radius=10)
    convex.propagate_ray(ray)
#    print(ray.p(),ray.k())
    output = g.OutputPlane(300)
    output.propagate_ray(ray)
    plt.xlabel('z (mm)')
    plt.ylabel('y (mm)')
    plt.grid()
    plt.plot(np.array(ray.vertices())[:, 2], np.array(ray.vertices())[:, 1])

# %%
# Tracing a few rays through the plano-convex lens reversed
for i in np.arange(-0.1, 0.12, 0.02):
    ray = r.Ray(pos=[0.0, i, 0.0], dirn=[0.0, 0.0, 1.0])
    convex = g.SphericalRefraction(z0=100, curvature=0.02, n1=1.0, n2=1.5168, aperture_radius=10)
    convex.propagate_ray(ray)
    plano = g.SphericalRefraction(z0=105, curvature=0.0, n1=1.5168, n2=1.0, aperture_radius=10)
    plano.propagate_ray(ray)
    output = g.OutputPlane(300)
    output.propagate_ray(ray)
    plt.xlabel('z (mm)')
    plt.ylabel('y (mm)')
    plt.grid()
    plt.plot(np.array(ray.vertices())[:, 2], np.array(ray.vertices())[:, 1])
