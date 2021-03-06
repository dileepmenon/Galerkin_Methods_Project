In this work, the linearized [shallow water equations](https://en.wikipedia.org/wiki/Shallow_water_equations) are discretized using the Continuous Galerkin (CG) method. All details 
about the discretization using linear triangular elements can be found at [Jose Abell's Research blog](http://www.joseabell.com/finite-elements-for-shallow-water-equations.html)

The rectangular domain is discretised into triangles using
[MeshPy](https://mathema.tician.de/software/meshpy/). The edge of the triangle
is chosen as 0.1 cm. The resulting mesh looks like ![Triangular mesh
](mesh2.png)

A sample simulation run for 0.165 sec with a gaussian perturbation applied at
the center of the tub is shown in the below animation. This work was animated using
using the 3D visualization tool [Mayavi](http://docs.enthought.com/mayavi/mayavi/)

![wave propogation](wave.gif)



