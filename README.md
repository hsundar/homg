# HoMG

High-order finite-element package using hexahedral elements. The code is a
testbed for geometric multigrid approaches for high order discretizations. The
current implementation supports setting up a combination of $h$ and $p$
heirarchy. The following smoothers are supported,
 * Jacobi
 * Chebyshev-accelerated Jacobi
 * block Jacobi
 * Symmetric SOR 

[Project Page](http://hsundar.github.io/homg/)

Details about the implementation and a comparison of the different methods can
be found in   

[Comparison of multigrid algorithms for high-order continuous finite element discretizations](http://onlinelibrary.wiley.com/doi/10.1002/nla.1979/abstract), _Hari Sundar, Georg Stadler, George Biros_, Numer. Linear Algebra Appl., doi: [10.1002/nla.1979](http://dx.doi.org/10.1002/nla.1979).

or
   
[Comparison of Multigrid Algorithms for High-order Continuous Finite Element
Discretizations](http://arxiv.org/pdf/1402.5938v1) _Hari Sundar, Georg Stadler, George Biros_
[arXiv](http://arxiv.org/abs/1402.5938)


## Basic Usage

A simple example in 2D
```matlab

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

```matlab

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

Using low-order preconditioning for high-order operator in CG method

```matlab

% specify the coefficients for a 2D example
mu = @(x,y)(1 + 10^6*( cos(2*pi*x)^2 + cos(2*pi*y)^2));

% create low-order approximation based on nodes of high-order operator
% we use a 2D unit square, discretized by 8 x 8 3rd-order elements
% usually algebraic MG is used to solve the low-order system
% here we use a direct solver

[it0, it1, it3] = low_order_precon([8 8], @homg.xform.identity, 3);

% Output is number of iterations with 0,1 or 3 Chebyshev smoothing steps
% on the finest level using the high-order residual

% And now for an example in 3D with warped geometry and 3x3x3 4th order elements

mu = @(x,y,z)(1 + 10^6*( cos(2*pi*x)^2 + cos(2*pi*y)^2 + cos(2*pi*z)^2) );
[it0, it1, it3] = low_order_precon([3 3 3], @homg.xform.shell, 4);


```
