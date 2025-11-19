# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 13:50:06 2021

@author: Gerry O'Brien
"""

import GeneralBaseClass as g
import RayClass as r
import numpy as np
import matplotlib.pyplot as plt


def parax_focus():  # Function finding the paraxial focal point
    """
    Function to find the paraxial focal point of the spherical lens
    """
    parax_ray = r.Ray(pos=[0, 0.1, 0], dirn=[0, 0, 1])
    refractor = g.SphericalRefraction()
    refractor.propagate_ray(parax_ray)
    l = -parax_ray.p()[1]/(parax_ray.k()[1]/np.linalg.norm(parax_ray.k()))
    z = parax_ray.p()[2]+l*(parax_ray.k()[2]/np.linalg.norm(parax_ray.k()))
    return z  # focal point


# Plotting a few rays through the spherical refractor
for i in np.arange(-0.1, 0.12, 0.02):
    ray = r.Ray(pos=[0.0, i, 0.0], dirn=[0.0, 0.0, 1.0])
    refractor = g.SphericalRefraction()
    refractor.propagate_ray(ray)
    output = g.OutputPlane(250)
    output.propagate_ray(ray)
    plt.title('A few rays through the spherical lens')
    plt.xlabel('z (mm)')
    plt.ylabel('y (mm)')
    plt.grid()
    plt.plot(np.array(ray.vertices())[:, 2], np.array(ray.vertices())[:, 1])
