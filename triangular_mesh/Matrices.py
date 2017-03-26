import numpy as np

from tri_mesh import mesh_rect


def tri_area(x_cord_list, y_cord_list):
    """Calculates the area of triangle"""
    arr = np.array([x_cord_list, y_cord_list, [1, 1, 1]])
    return 0.5*abs(np.linalg.det(arr))


def global_mass_matrix(mesh_points, mesh_elements):
    """Computes and returns the global mass matrix"""
    me_pt = mesh_points
    me_el = mesh_elements
    #Number of nodes in mesh
    num_nodes = len(mesh_points)
    global_matrix = np.zeros((num_nodes, num_nodes))
    for i in me_el:
        A = tri_area([me_pt[i[0]][0], me_pt[i[1]][0], me_pt[i[2]][0]],
                     [me_pt[i[0]][1], me_pt[i[1]][1], me_pt[i[2]][1]],
                    )
        #Mass matrix of element
        M_e = (A/12.0)*(np.mat('2 1 1; 1 2 1; 1 1 2'))
        for n, index in enumerate(i):
            row = global_matrix[index]
            #Other elements in list i excluding index
            oth_el = filter(lambda x: x!=index, i)
            a = M_e[n].item(n)
            row[index] += a
            row[oth_el[0]] += 0.5*a
            row[oth_el[1]] += 0.5*a
    return np.mat(global_matrix)


def global_K_hu_matrix(mesh_points, mesh_elements, H):
    """Computes the coefficient matrix of u in the weak form of the
    continuity equation

    Args:
        H : Depth of the water surface from undisturbed position
    """
    me_pt = mesh_points
    me_el = mesh_elements
    #Number of nodes in mesh
    num_nodes = len(mesh_points)
    global_matrix = np.zeros((num_nodes, num_nodes))
    for i in me_el:
        y_23, y_31, y_12 = (me_pt[i[1]][1] - me_pt[i[2]][1],
                            me_pt[i[2]][1] - me_pt[i[0]][1],
                            me_pt[i[0]][1] - me_pt[i[1]][1]
                           )
        m = np.matrix([[y_23, y_23, y_23],
                       [y_31, y_31, y_31],
                       [y_12, y_12, y_12]
                      ]
                     )
        #K_hu matrix of element
        K_hu_e = (H/6.0)*(m)
        for n, index in enumerate(i):
            #Other elements in list i excluding index
            oth_el = filter(lambda x: x!=index, i)
            row = global_matrix[index]
            a = K_hu_e[n].item(n)
            row[index] += a
            row[oth_el[0]] += a
            row[oth_el[1]] += a
    return -1*np.matrix(global_matrix)


def global_K_hv_matrix(mesh_points, mesh_elements, H):
    """Computes the coefficient matrix of v in the weak form of the
    continuity equation

    Args:
        H : Depth of the water surface from undisturbed position
    """
    me_pt = mesh_points
    me_el = mesh_elements
    #Number of nodes in mesh
    num_nodes = len(mesh_points)
    global_matrix = np.zeros((num_nodes, num_nodes))
    for i in me_el:
        x_32, x_13, x_21 = (me_pt[i[2]][0] - me_pt[i[1]][0],
                            me_pt[i[0]][0] - me_pt[i[2]][0],
                            me_pt[i[1]][0] - me_pt[i[0]][0]
                           )
        m = np.matrix([[x_32, x_32, x_32],
                       [x_13, x_13, x_13],
                       [x_21, x_21, x_21]
                      ]
                     )
        #K_hu matrix of element
        K_hu_e = (H/6.0)*(m)
        for n, index in enumerate(i):
            #Other elements in list i excluding index
            oth_el = filter(lambda x: x!=index, i)
            row = global_matrix[index]
            a = K_hu_e[n].item(n)
            row[index] += a
            row[oth_el[0]] += a
            row[oth_el[1]] += a
    return -1*np.matrix(global_matrix)
