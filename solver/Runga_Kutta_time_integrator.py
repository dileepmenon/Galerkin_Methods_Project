import math


def runge_kutta(dt, tf, M_inv, K_hu, K_hv, K_uh, K_vh, h_ini,
                u_ini, v_ini, index_LBP, list_index_CP,
                list_neigh_index):
    """Iterates the continuity and momentum equations in time to find the
    updated values of [h, u, v]

    Args:
        dt                : Timestep
        tf                : Final time to which simulation is run
        M_inv             : Inverse of mass matrix
        K_hu              : Coefficient matrix of u in continuity equation
        K_hv              : Coefficient matrix of v in continuity equation
        K_uh              : Coefficient matrix of h in x-momentum equation
        K_vh              : Coefficient matrix of h in y-momentum equation
        h_ini             : Initial vector of h values in nodes
        u_ini             : Initial vector of u values in nodes
        v_ini             : Initial vector of v values in nodes
        index_LBP         : The index of the last boundary node
        list_index_CP     : List containing the index of corner points
        list_neigh_index  : List containing the indices of the nearest
                            node of the boundary nodes

    Returns:
        h_list            : List of updated values of h in each iteration
        u_list            : List of updated values of u in each iteration
        v_list            : List of updated values of v in each iteration
    """
    #Number of iterations
    n = int(tf/dt)

    #Coefficient matrices
    A = -1*M_inv*K_uh
    B = -1*M_inv*K_vh
    P = M_inv*K_hu
    Q = M_inv*K_hv

    h_list = []
    u_list = []
    v_list = []

    #Initial values
    h = h_ini
    u = u_ini
    v = v_ini

    #Append initial values of [h, u, v] to respective lists
    h_list.append(h)
    u_list.append(u)
    v_list.append(v)

    #Note: Boundary Condition Implementation
    #The Boundary condition for reflection is du/dn = 0 for left and right
    #boundary, and dv/dn = 0 for top and bottom boundary, where n is normal
    #vector. To ensure this, all boundary points except the corner points
    #should have u/v value equal to u/v of nearest neighbour perpendicular to
    #their edge. For corner points it should be same as diagonal point.

    #Starting and ending index of all boundary points excluding corner points
    bt_edge_index = [list_index_CP[0]+1, list_index_CP[1]-1] #Bottom edge
    rt_edge_index = [list_index_CP[1]+1, list_index_CP[2]-1] #Right edge
    tp_edge_index = [list_index_CP[2]+1, list_index_CP[3]-1] #Top edge
    lt_edge_index = [list_index_CP[3]+1, index_LBP] #Left edge

    #Indices of the nearest nodes of boundary points of respective edge
    rt_nei = list_neigh_index[rt_edge_index[0]:rt_edge_index[1]+1] #Right edge
    lt_nei = list_neigh_index[lt_edge_index[0]:lt_edge_index[1]+1] #Left edge
    tp_nei = list_neigh_index[tp_edge_index[0]:tp_edge_index[1]+1] #Top edge
    bt_nei = list_neigh_index[bt_edge_index[0]:bt_edge_index[1]+1] #Bottom edge

    print n
    for i in range(n):
        print i
        #Values of [h, u, v] after half time step
        h_n_half = h + (dt/2.0)*((P*u)+(Q*v))
        u_n_half = u + (dt/2.0)*A*h
        v_n_half = v + (dt/2.0)*B*h

        #To make sure du/dx = 0 at left boundary
        u_n_half_lt_nei = [u_n_half[i] for i in lt_nei]
        u_n_half[lt_edge_index[0]:lt_edge_index[1]+1] = u_n_half_lt_nei[:]

        #To make sure du/dx = 0 at right boundary
        u_n_half_rt_nei = [u_n_half[i] for i in rt_nei]
        u_n_half[rt_edge_index[0]:rt_edge_index[1]+1] = u_n_half_rt_nei[:]

        #To make sure dv/dy = 0 at top boundary
        v_n_half_tp_nei = [v_n_half[i] for i in tp_nei]
        v_n_half[tp_edge_index[0]:tp_edge_index[1]+1] = v_n_half_tp_nei[:]

        #To make sure dv/dy = 0 at bottom boundary
        v_n_half_bt_nei = [v_n_half[i] for i in bt_nei]
        v_n_half[bt_edge_index[0]:bt_edge_index[1]+1] = v_n_half_bt_nei[:]

        #To make sure du/dx and dv/dy = 0 at corner points
        ind_cnr_pt_1 = list_index_CP[0]
        u_n_half[ind_cnr_pt_1] = u_n_half[list_neigh_index[ind_cnr_pt_1]]
        v_n_half[ind_cnr_pt_1] = v_n_half[list_neigh_index[ind_cnr_pt_1]]

        ind_cnr_pt_2 = list_index_CP[1]
        u_n_half[ind_cnr_pt_2] = u_n_half[list_neigh_index[ind_cnr_pt_2]]
        v_n_half[ind_cnr_pt_2] = v_n_half[list_neigh_index[ind_cnr_pt_2]]

        ind_cnr_pt_3 = list_index_CP[2]
        u_n_half[ind_cnr_pt_3] = u_n_half[list_neigh_index[ind_cnr_pt_3]]
        v_n_half[ind_cnr_pt_3] = v_n_half[list_neigh_index[ind_cnr_pt_3]]

        ind_cnr_pt_4 = list_index_CP[3]
        u_n_half[ind_cnr_pt_4] = u_n_half[list_neigh_index[ind_cnr_pt_4]]
        v_n_half[ind_cnr_pt_4] = v_n_half[list_neigh_index[ind_cnr_pt_4]]

        #Values of [h, u, v] after one timestep
        h_n = h + dt*((P*u_n_half)+(Q*v_n_half))
        u_n = u + dt*A*h_n_half
        v_n = v + dt*B*h_n_half

        #To make sure du/dx = 0 at left boundary
        u_n_lt_nei = [u_n[i] for i in lt_nei]
        u_n[lt_edge_index[0]:lt_edge_index[1]+1] = u_n_lt_nei[:]

        #To make sure du/dx = 0 at right boundary
        u_n_rt_nei = [u_n[i] for i in rt_nei]
        u_n[rt_edge_index[0]:rt_edge_index[1]+1] = u_n_rt_nei[:]

        #To make sure dv/dy = 0 at top boundary
        v_n_tp_nei = [v_n[i] for i in tp_nei]
        v_n[tp_edge_index[0]:tp_edge_index[1]+1] = v_n_tp_nei[:]

        #To make sure dv/dy = 0 at bottom boundary
        v_n_bt_nei = [v_n[i] for i in bt_nei]
        v_n[bt_edge_index[0]:bt_edge_index[1]+1] = v_n_bt_nei[:]

        #To make sure du/dx and dv/dy = 0 at corner points
        u_n[ind_cnr_pt_1] = u_n[list_neigh_index[ind_cnr_pt_1]]
        v_n[ind_cnr_pt_1] = v_n[list_neigh_index[ind_cnr_pt_1]]

        u_n[ind_cnr_pt_2] = u_n[list_neigh_index[ind_cnr_pt_2]]
        v_n[ind_cnr_pt_2] = v_n[list_neigh_index[ind_cnr_pt_2]]

        u_n[ind_cnr_pt_3] = u_n[list_neigh_index[ind_cnr_pt_3]]
        v_n[ind_cnr_pt_3] = v_n[list_neigh_index[ind_cnr_pt_3]]

        u_n[ind_cnr_pt_4] = u_n[list_neigh_index[ind_cnr_pt_4]]
        v_n[ind_cnr_pt_4] = v_n[list_neigh_index[ind_cnr_pt_4]]

        #Update values of [h, u, v]
        h = h_n
        u = u_n
        v = v_n

        #Append updated values after each timestep to respective list
        h_list.append(h)
        u_list.append(u)
        v_list.append(v)

    return h_list, u_list, v_list
