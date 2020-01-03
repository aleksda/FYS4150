from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
from numba import njit, prange
import numpy as np

"""
This source code contains two classes for simulating a solar system. 
The first is called Planet and it contains information on the planet of choise. 

The second is called System and is used for solving and plotting the equations needed for the system.
The System class contains two functions: solve and plot. 
The solve function has four inner functions: a_g, forward_euler, euler_cromer and velocity_verlet.
The a_g calculates the acceleration and is used by the forward_euler, euler_cromer and velocity_verlet functions to calculate the next pos.
"""

class Planet:

    def __init__(self, x0, v0, m, name):

        self.x0 = np.array(x0)
        self.v0 = np.array(v0)
        self.m = m
        self.name = name

        self.G = 6.673e-11

        # For earth/sun
        self.MASS   = 1989.1
        self.RADIUS = 1.4766957e+12

    def set(self, t, x, v):

        self.t = t
        self.x = x
        self.v = v

    def get_escape_velocity(self):

        return np.sqrt((2 * self.G * self.MASS) / self.RADIUS)

Mercury = Planet(x0 = [-0.06212132,  0.31238972,  0.03030863],
                    v0 = [-12.14754017,   1.5667267,    0.98618233],
                    m = 1.65919050e-07, name = "Mercury")

Venus = Planet(x0 = [ 0.33347755,  0.63727269, -0.02827693],
                v0 = [ 6.49052304,  3.3935057,  -0.32807182],
                m = 2.44856295e-06, name = "Venus")

Earth = Planet(x0 = [ 5.33585891E-1,  8.37073194E-1, -2.81264913E-5],
                v0 = [-5.37578469,  3.39044094, -2.10466761e-04],
                m = 3.00162645e-06, name = "Earth")

"""
Earth = Planet(x0 = [ 5.33585891E-1,  8.37073194E-1, -2.81264913E-5],
                v0 = [-5.37578469,  2*np.sqrt(2)*np.pi, -2.10466761e-04],
                m = 3.00162645e-06, name = "Earth")                
"""
"""
Earth = Planet(x0 = [9.413801075750535E-01, 3.379019986046322E-01, 
                        -9.334104672733438E-05],
                v0 = [-5.994522787486753E-03, 1.617377250092178E-02,
                        -1.732657683299539E-07],
                m = 5.97219e24, name = "Earth")
"""
Mars = Planet(x0 = [-1.58321451,  0.39412134,  0.03035798],
                v0 = [ 1.4458615,   4.51404483, -0.13004533],
                m = 3.22787970e-07, name = "Mars")

Jupiter = Planet(x0 = [0.21011067, 5.23090799, 0.01699342],
                    v0 = [ 2.7184775,   0.24193614, -0.06181532],
                    m = 9.54285930e-04, name = "Jupiter")
"""
Jupiter = Planet(v0 = [-2.666952709077877E+00, -4.655671225645230E+00,
                        7.896515774211305E-02],
                    x0 = [6.458958874387921E-03, -3.390642961368397E-03,
                        -1.303431975919576E-04],
                    m = 1898.13e24, name = "Jupiter")
"""
Saturn = Planet(x0 = [3.58856249, 9.36621316, 0.02000048],
                v0 = [ 1.7885416,   0.72238731, -0.08369087],
                m = 2.85581880e-04, name = "Saturn")

Uranus = Planet(x0 = [16.31727197, 11.25858689, -0.16957772],
                v0 = [-0.82585914,  1.11470731,  0.01480277],
                m = 4.36417380e-05, name = "Uranus")

Neptune = Planet(x0 = [29.21158221,  6.48930242, -0.53957529],
                    v0 = [ 0.24095139,  1.12523358, -0.02887421],
                    m = 5.12840700e-05, name = "Neptune")

Pluto = Planet(x0 = [12.84801582, 31.38483364, -0.35803234],
                v0 = [ 1.08913103,  0.19086789, -0.33884746],
                m = 7.34066100e-09, name = "Pluto")

print(Earth.get_escape_velocity())