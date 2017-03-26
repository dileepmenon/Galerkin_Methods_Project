import math


def runge_kutta(dt, tf, M_inv, K_hu, K_hv, K_uh, K_vh, h_ini, u_ini, v_ini):
    """Iterates the continuity and momentum equations in time to find the
    updated values of [h, u, v]

    Args:
        dt      : Timestep
        tf      : Final time to which simulation is run
        M_inv   : Inverse of mass matrix
        K_hu    : Coefficient matrix of u in continuity equation
        K_hv    : Coefficient matrix of v in continuity equation
        K_uh    : Coefficient matrix of h in x-momentum equation
        K_vh    : Coefficient matrix of h in y-momentum equation
        h_ini   : initial vector of h values in nodes
        u_ini   : initial vector of u values in nodes
        v_ini   : initial vector of v values in nodes

    Returns:
        h_list  : list of updated values of h in each iteration
        u_list  : list of updated values of u in each iteration
        v_list  : list of updated values of v in each iteration
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

    for i in range(n):
        #Values of [h, u, v] after half time step
        h_n_half = h + (dt/2.0)*((P*u)+(Q*v))
        u_n_half = u + (dt/2.0)*A*h
        v_n_half = v + (dt/2.0)*B*h

        #Values of [h, u,v ] after one timestep
        h_n = h + dt*((P*u_n_half)+(Q*v_n_half))
        u_n = u + dt*A*h_n_half
        v_n = v + dt*B*h_n_half

        #Update values of [h, u, v]
        h = h_n
        u = u_n
        v = v_n

        #Append updated values after each timestep to respective list
        h_list.append(h)
        u_list.append(u)
        v_list.append(v)

    return h_list, u_list, v_list
