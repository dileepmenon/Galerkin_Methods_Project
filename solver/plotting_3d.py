from mayavi import mlab

from solver import solve

from tri_mesh import mesh_rect

def plot_and_animate():
    h, u, v = solve()
    mp, me, l = mesh_rect()
    x = np.array([i[0] for i in mp])
    y = np.array([i[1] for i in mp])
    src = mlab.pipline.scalar_scatter(x, y, np.array(h[0].T)[0],
                                      np.array(h[0].T)[0]
                                     )
    f = mlab.pipeline.delaunay2d(src)
    s = mlab.pipeline.surface(f)
    mlab.outline()
    mlab.axes(src)
    mlab.colorbar()
    for i in range(1000):
        s.mlab_source.set(z = np.array(h[i].T)[0], scalars = np.array(h[i].T)[0])
        mlab.savefig('/tmp/img%04d.png'%i)

if __name__ == '__main__':
    plot_and_animate()

