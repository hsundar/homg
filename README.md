# HoMG

High-order finite-element package using hexahedral elements. The code is a
testbed for geometric multigrid approaches for high order discretizations. The
current implementation supports setting up a combination of $h$ and $p$
heirarchy. The following smoothers are supported,
 * Jacobi
 * Chebyshev-accelerated Jacobi
 * block Jacobi
 * Symmetric SOR 


## Basic Usage

A simple example in 2D
```
#!matlab

% specify the coefficients
mu = @(x,y)(1 + 1e6*( cos(2*pi*x)^2 + cos(2*pi*y)^2 ) );

% create the mesh heirarchy, for a warped mesh
% in this case we create a h+p heirarchy
%     grid 4 --> 32x32, p=4  (finest)
%     grid 3 --> 32x32, p=2      
%     grid 2 --> 32x32, p=1
%     grid 1 --> 16x16, p=1
%     grid 0 -->  8x8,  p=1  (coarsest)

g = create_hexmesh_grids(2, mu, @homg.xform.shell, [1 2 4], [8 16 32]);

% now solve using multigrid as the solver and the choice of smoother
g.solve (150, 'jacobi', 3,3, g.L, g.get_u0);
g.solve (150, 'chebyshev', 3,3, g.L, g.get_u0);
g.solve (150, 'ssor', 2,1, g.L, g.get_u0);

% or solve using CG preconditioned using multigrid
g.solve_pcg(150, 'jacobi', 3,3, g.L, g.get_u0);
g.solve_pcg(150, 'chebyshev', 3,3, g.L, g.get_u0);
g.solve_pcg(150, 'ssor', 2,1, g.L, g.get_u0);

```
The 3D example is similar with a few changes in the grid setup.

```
#!matlab

% specify the coefficients
mu = @(x,y,z)(1 + 10^6*( cos(2*pi*x)^2 + cos(2*pi*y)^2 + cos(2*pi*z)^2) );

% create the mesh heirarchy, for a warped mesh
% in this case we create a h+p heirarchy
%     grid 4 --> 8x8x8, p=4  (finest)
%     grid 3 --> 8x8x8, p=2      
%     grid 2 --> 8x8x8, p=1
%     grid 1 --> 4x4x4, p=1
%     grid 0 --> 2x2x2, p=1  (coarsest)

g = create_hexmesh_grids(3, mu, @homg.xform.shell, [1 2 4], [2 4 8]);

% now solve using multigrid as the solver and the choice of smoother
g.solve (150, 'jacobi', 3,3, g.L, g.get_u0);
g.solve (150, 'chebyshev', 3,3, g.L, g.get_u0);
g.solve (150, 'ssor', 2,1, g.L, g.get_u0);

% or solve using CG preconditioned using multigrid
g.solve_pcg(150, 'jacobi', 3,3, g.L, g.get_u0);
g.solve_pcg(150, 'chebyshev', 3,3, g.L, g.get_u0);
g.solve_pcg(150, 'ssor', 2,1, g.L, g.get_u0);

```

