# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 16:23:39 2021

@author: Gerry O'Brien
"""
import numpy as np


class Ray:

    """
    A class that represents an optical ray
    It can be used to generate a series of points the ray path follows
    """

    def __init__(self, pos=[0.0, 0.0, 0.0], dirn=[0.0, 0.0, 0.0]):  # initialises at 0,0,0 with 1,0,0 dirn
        self._p = np.array(pos)
        self._k = np.array(dirn)
        self._v = [pos]

    def p(self):
        return self._p

    def k(self):
        return self._k

    def append(self, p=[0.0, 0.0, 0.0], k=[0.0, 0.0, 0.0]):  # adds a new point with a new direction and add pos to vertices list
        self._p = np.array(p)
        self._k = np.array(k)
        self._v.append(p)

    def vertices(self):
        return self._v

# TESTS ON INITIALISING, FUNCTIONS WORK, PROPERLY RECORDS POINTS
# ray=Ray(pos=[1,2,3],dirn=[4,5,6])
# print(ray.p())
# print(ray.k())
# ray.append([2,2,2],[1,1,1])
# ray.append([5,6,7],[9,8,7])
# print(ray.vertices())
