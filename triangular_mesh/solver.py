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
        h_ini.append(0.1*(1/(2*np.pi*sigma_x*sigma_y) * np.exp(-(x**2/(2*sigma_x**2)
                     + y**2/(2*sigma_y**2)))))

    return np.mat(h_ini).T


def least_dist_point_index(mesh_point_i, list_mesh_tup):
    """Finds the index of the nearest node from node i by comparing
       the distances of the neighbouring nodes.

    Args:
        mesh_point_i  : (x,y) co-ordinates of the node i
        list_mesh_tup : list of tuples containing ((x,y), index) of
                        neighbouring nodes
    """
    #Stores the ((x-x_i)**2+(y-y_i)**2, index) of the neighbouring nodes of
    #node i
    dist = []

    for tup in list_mesh_tup:
        dist.append(((tup[0][0] - mesh_point_i[0])**2 \
                     + (tup[0][1] - mesh_point_i[1])**2,
                     tup[1]
                    )
                   )

    #Nodes are sorted in increasing distance from node i, the first element
    #contains the ((x-x_i)**2+(y-y_i)**2, index) of the nearest node
    index_nearest_neighbour = sorted(dist, key= lambda d : d[0])[0][1]

    return index_nearest_neighbour


def nearest_neighbour(index_LBP, mesh_points, mesh_elements):
    """Finds the index of the nearest node of all boundary nodes.

    Args:
        index_LBP      : The index of the last boundary node
        mesh_points    : List containing the (x,y) co-ordinates of all mesh
                         nodes
        mesh_elements  : List containing the node indices (n1, n2, n3) of all
                         triangles in the mesh

    Returns:
        list_neighbour : list containing the index of all the nearest nodes of
                         the boundary nodes

    """
    list_neighbour = []
    for ele_num in range(index_LBP+1):
        #Filters the triangluar elements which contains node ele_num
        list_tri = filter(lambda element: element[0] == ele_num
                          or element[1] == ele_num
                          or element[2] == ele_num,
                          mesh_elements
                         )
        #From the triangle elements containing node ele_num, stores the
        #unique indices of the nodes (Note: Also contains node ele_num)
        list_node = set([node_num for tri in list_tri for node_num in tri])

        #List of tuples containing ((x,y), index) of neighbouring nodes
        list_mesh_tup = []

        for node_num in list_node:
            if (node_num != ele_num
                and mesh_points[node_num][0] != mesh_points[ele_num][0]
                and mesh_points[node_num][1] != mesh_points[ele_num][1]
               ):

                list_mesh_tup.append((mesh_points[node_num], node_num))

        list_neighbour.append(least_dist_point_index(mesh_points[ele_num],
                                                     list_mesh_tup
                                                    )
                             )

    return list_neighbour[:]


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
    tf = 1.0

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
    u_ini = np.mat(np.zeros((len(mesh_points), 1)))
    v_ini = np.mat(np.zeros((len(mesh_points), 1)))

    #To find the boundary points
    #In the list mesh_points, the boundary points are stored first and then the
    #interior points. The boundary points in the bottom edge of the domain
    #comes first and then in a anticlockwise fashion the right, top and left
    #edge points are stored. To find index of last boundary boundary point we
    #make use of the fact that the value of (x,y) of that point is zero and non-
    #zero respectively and value of x of the next point should be non zero.
    for i, r in enumerate(mesh_points):
        if r[0] == 0 and r[1] != 0 and mesh_points[i+1][0] != 0:
            index_LBP = i # Index of last boundary point

    # Stores the position [x, y] of the corner points in domain
    pos_CP = [[0., 0.], [5., 0.], [5., 2.5], [0., 2.5]]

    arr_eq = np.array_equal

    # List which stores the index of the corner points of mesh_points array
    list_index_CP = [ind for ind, r in enumerate(mesh_points)
                     if arr_eq(r, pos_CP[0]) or arr_eq(r, pos_CP[1])
                     or arr_eq(r, pos_CP[2]) or arr_eq(r, pos_CP[3])
                    ]

    #list containing the (x,y) values of all boundary points
    #boun_points = mesh_points[:index_LBP+1]

    #list containing the index of nearest neighbor of each boundary point
    list_neigh_index = nearest_neighbour(index_LBP, mesh_points, mesh_elements)

    #Time Integration
    h_list, u_list, v_list = runge_kutta(dt, tf, M_inv,
                                         K_hu, K_hv,
                                         K_uh, K_vh,
                                         h_ini, u_ini,
                                         v_ini, index_LBP,
                                         list_index_CP,
                                         list_neigh_index
                                        )

    return h_list, u_list, v_list
