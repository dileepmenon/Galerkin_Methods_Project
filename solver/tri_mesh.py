import jw_meshtools as mt


def mesh_rect(tri_length = 0.1, x_dim = 5.0, y_dim = 2.5):
    """ Creates a triangular mesh on a rectangular domain with co-ordinates :
    [(0, 0), (5, 0), (5, 2.5), (0, 2.5)]
    """
    #length of edge of triangle
    length = tri_length
    p, v = mt.RectangleSegments([0, 0], [x_dim, y_dim], edge_length = length)
    mesh_points, mesh_elements = mt.DoTriMesh(p, v, edge_length = length)
    return mesh_points, mesh_elements, length


