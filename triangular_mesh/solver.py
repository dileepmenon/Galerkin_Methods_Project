import numpy as np

import matplotlib.pyplot as plt

import math

from Matrices import *

from tri_mesh import mesh_rect

from Runga_Kutta_time_integrator import runge_kutta


def gaussian(mesh_points):
    """Gives a gaussian function to get the initial height distribution

    Args:
        x : list of x co-ordinate of mesh
        y : list of y co-ordinate of mesh

    Returns:
        Column vector of h values at each node of mesh
    """
    #Standard deviations
    sigma_x = 0.2
    sigma_y = 0.2

    h_ini = []

    for x,y in mesh_points:
        #Transform co-ordinate to get gaussian in center of rectangular domain
        x = x - 2.5
        y = y - 1.25

        #Height distribution at each node as a gaussian
        h_ini.append(1/(2*np.pi*sigma_x*sigma_y) * np.exp(-(x**2/(2*sigma_x**2)
                     + y**2/(2*sigma_y**2))))

    return np.mat(h_ini).T


def solve():
    g = 9.8 #acceleration due to gravity in m/s^2
    H = 20 #Height of undistrubed surface in cm

    #Triangle mesh points, element indices and edge length in cm
    mesh_points, mesh_elements, dx = mesh_rect()

    #Wave speed
    c = math.sqrt(g*(H/100.0))

    #Courant Number((c*dt)/dx) for stability of runge kutta should be <= 1/4
    courant_num = 1/4.0

    #Time Step
    dt = (dx*0.01*courant_num)/c

    #Final time in sec
    tf = 5.0

    #Mass Matrix
    M = global_mass_matrix(mesh_points, mesh_elements)

    #Inverse of Mass Matrix
    M_inv = M.I

    #Global K_hu matrix
    K_hu =  global_K_hu_matrix(mesh_points, mesh_elements, H)

    #Global K_hv matrix
    K_hv =  global_K_hv_matrix(mesh_points, mesh_elements, H)

    #Global K_uh matrix
    K_uh = global_K_uh_matrix(mesh_points, mesh_elements)

    #Global K_vh matrix
    K_vh = global_K_vh_matrix(mesh_points, mesh_elements)

    #Gaussian Initial Condition for height
    h_ini = gaussian(mesh_points)

    #Initial condition for velocities in x and y direction
    u_ini = np.mat(np.zeros((240,1)))
    v_ini = np.mat(np.zeros((240,1)))

    #Boundary Conditions

    #Time Integration
    h_list, u_list, v_list = runge_kutta(dt, tf, M_inv,
                                         K_hu, K_hv,
                                         K_uh, K_vh,
                                         h_ini, u_ini,
                                         v_ini
                                        )

    return h_list, u_list, v_list
