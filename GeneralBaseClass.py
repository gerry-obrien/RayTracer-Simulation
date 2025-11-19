# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 15:30:28 2021

@author: Gerry O'Brien
"""
import numpy as np
import matplotlib.pyplot as plt


class OpticalElement:
    """
    Class to propagate a ray
    """

    def propagate_ray(self, ray):
        raise NotImplementedError()


class SphericalRefraction(OpticalElement):
    """
    An inherited class to model a spherical lens
    Parameters:
        z0 = the position of the surface along the z-axis
        curvature = 1/radius of sphere
        n1 = refractive index before the boundary
        n2 = refractive index after the boundary
    """

    def __init__(self, z0=100, curvature=0.03, n1=1.0, n2=1.5, aperture_radius=10):
        self._z0 = z0
        self._curv = curvature
        self._n1 = n1
        self._n2 = n2
        self._ar = aperture_radius

    def intercept(self, ray):
        """
        A function to find the intercept with the spherical surface
        """
        k = ray.k()
        k_hat = k/np.sqrt(k[0]**2+k[1]**2+k[2]**2)
        # k is the direction vector from start point to optical aparatus
        if self._curv == 0.0:
            l = (self._z0 - ray.p()[2])/(ray.k()[2]/np.sqrt(k[0]**2+k[1]**2+k[2]**2))
            intercept = ray.p() + l*k_hat
            return intercept

        R = 1/self._curv  # radius
        O_z = self._z0 + R  # centre of sphere
        r = ray.p() - np.array([0, 0, O_z])
        mod_r = np.sqrt(r[0]**2+r[1]**2+r[2]**2)
        # r is the vector from the centre of the sphere
        # to the initial ray position
        if (mod_r**2-R**2) < (np.dot(r, k_hat))**2:  # if there's intersection
            l1 = -np.dot(r, k_hat)+np.sqrt((np.dot(r, k_hat))**2-(mod_r**2-R**2))
            l2 = -np.dot(r, k_hat)-np.sqrt((np.dot(r, k_hat))**2-(mod_r**2-R**2))
            if l1 < 0 or l2 < 0:
                intercept = max(l1, l2) * k_hat + ray.p()
            if l1 > 0 and l2 > 0:
                intercept = min(l1, l2) * k_hat + ray.p()
            return intercept
        else:
            return None

    def SnellsLaw(self, ray):  # returns direction unit vector of refracted ray
        """
        Snell's Law finds the direction vector after the refraction
        """
        if self._curv == 0.0:
            # for a flat surface
            n = np.array([0, 0, 1])
            mu = self._n1/self._n2
            # k1 is the direction before refraction
            k1 = ray.k()/np.sqrt(ray.k()[0]**2+ray.k()[1]**2+ray.k()[2]**2)
            k2 = mu*k1 + n*(np.sqrt(1-mu**2*(1-np.dot(n, k1)**2))) - mu*n*np.dot(n, k1)
            return k2
        if self._curv < 0.0:
            # for negative curve
            R = 1/self._curv
            O_z = self._z0 + R
            mu = self._n1/self._n2
            k1 = ray.k()
            k1_hat = k1/np.sqrt(k1[0]**2+k1[1]**2+k1[2]**2)
            # n is the normal to the refractor surface
            n = - np.array([0, 0, O_z]) + self.intercept(ray)
            n_hat = n/np.sqrt(n[0]**2+n[1]**2+n[2]**2)
            k2_hat = mu*k1_hat + n_hat*(np.sqrt(1-mu**2*(1-np.dot(n_hat, k1_hat)**2))) - mu*n_hat*np.dot(n_hat, k1_hat)
            sin_theta = np.sin(np.arccos(np.dot(k1_hat, n_hat)))
            if sin_theta > 1/mu:  # conditions for total internal reflection
                return None
            else:
                return k2_hat
        if self._curv > 0.0:
            # for positive curve
            R = 1/self._curv
            O_z = self._z0 + R
            mu = self._n1/self._n2
            k1 = ray.k()
            k1_hat = k1/np.sqrt(k1[0]**2+k1[1]**2+k1[2]**2)
            # n is the normal to the refractor surface
            n = np.array([0, 0, O_z]) - self.intercept(ray)
            n_hat = n/np.sqrt(n[0]**2+n[1]**2+n[2]**2)
            k2_hat = mu*k1_hat + n_hat*(np.sqrt(1-mu**2*(1-np.dot(n_hat, k1_hat)**2))) - mu*n_hat*np.dot(n_hat, k1_hat)
            sin_theta = np.sin(np.arccos(np.dot(k1_hat, n_hat)))
            if sin_theta > 1/mu:  # conditions for total internal reflection
                return None
            else:
                return k2_hat

#    def output_plane(self, ray):
        # original output plane was a method within a class
#        z = self._z0 + self._n2/(self._curv*(self._n2-self._n1))
#        # z = 500
#        k = self.SnellsLaw(ray)
#        l = (z-self.intercept(ray)[2])/(k[2]/np.sqrt(k[0]**2+k[1]**2+k[2]**2))
#        k_hat = k/np.sqrt(k[0]**2+k[1]**2+k[2]**2)
#        output_intercept = self.intercept(ray) + l*k_hat
#        return output_intercept

    def propagate_ray(self, ray):
        """
        This propagates the ray after passing through the refractor
        """
        if type(self.intercept(ray)) == type(None):
            raise TypeError("There's no intercept")
#        print(self.SnellsLaw(ray))
        if type(self.SnellsLaw(ray)) == type(None):
            raise TypeError("No refraction as total internal reflection occurs")
        else:
            new_p = self.intercept(ray)
            # print(new_p)
            new_k = self.SnellsLaw(ray)
#            output_p = self.output_plane(ray)
            ray.append(new_p.tolist(), new_k)
            # had to change new position to a list to plot the vertices
#            ray.append(output_p.tolist(), new_k)
#        return ray


class OutputPlane(OpticalElement):
    """
    A class to find the points at and propagate to a specified output plane
    """

    def __init__(self, z):
        self._z = z

    def intercept(self, ray):
        """
        Method to find the intercept of a ray with the ouput plane
        """
        l = (self._z-ray.p()[2])/(ray.k()[2]/np.linalg.norm(ray.k()))
        output_intercept = ray.p()+l*(ray.k()/np.linalg.norm(ray.k()))
        return output_intercept

    def propagate_ray(self, ray):
        """
        Method to propagate the ray to the output plane
        """
        ray.append(self.intercept(ray).tolist(), ray.k().tolist())


def bundlefunc(R, n_step):
    """
    Function to create a bundle of rays with two arguments:
    R is the radius of the bundle,
    n_step is the number of concentric circles you want within the bundle
    """
    r_step = R/n_step
    theta_step = np.pi/n_step
    coord_list = []
    for i in np.arange(1, n_step+1):
        r = i*r_step
        for j in np.arange(1, 2*n_step+1):
            theta = j*theta_step
            x = r*np.cos(theta)
            y = r*np.sin(theta)
            z = 0
            coord_list.append([x, y, z])
    return coord_list
