Simulation example (Figure 1 in arXiv:2510.08443):
- Choose domain by commenting out the macro in the top of 2d-sfem-example/main.cpp.
- Set model coefficients in 2d-sfem-example/example-1.h (default is example from Figure 1).
- Set the value of gamma in 2d-sfem-example/main.cpp.
- Optional in 2d-sfem-example/main.cpp: set the parameter Nt (default 10) and NN_reference (default 6). Number of time steps i 2^Nt and the mesh is iteratively refined NN_reference times so that h \leq C 4^(-NN_reference).
- Run 2d-sfem-example/main.cpp.
- The result is saved in figures/motivation-transformed-sphere-(...).vtk where (...) depends on the choice of gamma and NN_reference.
- Open the .vtk file in Paraview.

1d convergence rate example (first row in Figure 2 in arXiv:2510.08443)
- Set the domain by defining a macro from 1d-sfem-rates/2dmesh.h in 1d-sfem-rates/main.cpp (e.g. write #define sphere_mesh at the top of 1d-sfem-rates/main.cpp, the default is line_mesh).
- Change b_process to 1. in 1d-sfem-rates/compute-error.h (this si from arXiv:2508.16458).
- Change gamma in 1d-sfem-rates/main.cpp.
- Run 1d-sfem-rates/main.cpp.

2d convergence rate example (second row in Figure 2 in arXiv:2510.08443)
- Set the domain by defining a macro from 2d-sfem-rates/2dmesh.h in 2d-sfem-rates/main.cpp (e.g. write #define icosasphere_mesh at the top of 2d-sfem-rates/main.cpp, the default is square_mesh).
- Change b_process to 1. in 1d-sfem-rates/compute-error.h (this si from arXiv:2508.16458).
- Change gamma in 2d-sfem-rates/main.cpp.
- Run 2d-sfem-rates/main.cpp.

1d convergence rate example (first row in Figure 1 in arXiv:2508.16458)
- Change gamma in 1d-sfem-rates/main.cpp.
- Run 1d-sfem-rates/main.cpp.

2d convergence rate example (second row in Figure 1 in arXiv:2508.16458)
- Change gamma in 2d-sfem-rates/main.cpp.
- Run 2d-sfem-rates/main.cpp.
