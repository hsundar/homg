% specify the coefficients
mu = @(x,y)(1 + 1e6*( cos(2*pi*x)^2 + cos(2*pi*y)^2 ) );

% create the mesh heirarchy, for a warped mesh
% in this case we create a h+p heirarchy
%     grid 4 --> 16x16, p=4  (finest)
%     grid 3 --> 16x16, p=2      
%     grid 2 --> 16x16, p=1
%     grid 1 -->  8x8 , p=1  (coarsest)

disp('--------------------');
disp('Creating Mesh heirarchy');
g = create_hexmesh_grids(2, mu, @homg.xform.shell, [1 2 4], [8 16]);

disp('--------------------');
disp('==== Solving using Multigrid as a solver ====');
% now solve using multigrid as the solver and the choice of smoother
disp('==== Jacobi ====');
g.solve (150, 'jacobi', 3,3, g.L, g.get_u0);
disp('==== Chebyshev ====');
g.solve (150, 'chebyshev', 3,3, g.L, g.get_u0);
disp('==== SSOR ====');
g.solve (150, 'ssor', 2,1, g.L, g.get_u0);

disp('--------------------');
disp('==== Solving using Multigrid as a preconditioner ====');
% or solve using CG preconditioned using multigrid
disp('==== Jacobi ====');
g.solve_pcg(150, 'jacobi', 3,3, g.L, g.get_u0);
disp('==== Chebyshev ====');
g.solve_pcg(150, 'chebyshev', 3,3, g.L, g.get_u0);
disp('==== SSOR ====');
g.solve_pcg(150, 'ssor', 2,1, g.L, g.get_u0);
disp('--------------------');

