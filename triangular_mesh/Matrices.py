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
            row[index], row[oth_el[0]], row[oth_el[1]] = (a,
                                                          0.5*a,
                                                          0.5*a
                                                          )
    return np.mat(global_matrix)

