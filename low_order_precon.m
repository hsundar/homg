% Example usage: low_order_precon([16 16], @homg.xform.identity, 3)
function [it_gll, it_smooth_1, it_smooth_3] = low_order_precon(nelems, xform, order)

mesh = homg.hexmesh(nelems, xform); 

if (mesh.dim == 2)
  mu = @(xx,yy)(1 + 1000000*( (cos(2*pi*xx))^2 + (cos(2*pi*yy))^2 ));
else
  mu = @(xx,yy,zz)(1 + 1000000*( (cos(2*pi*xx))^2 + (cos(2*pi*yy))^2 + (cos(2*pi*zz))^2 ));
end

mesh.set_coeff(mu);

[K, M]          =  mesh.assemble_poisson (order);
[K_lin, M_lin]  =  mesh.assemble_poisson_linearized (order);

grid_lin = homg.grid(mesh, order);
grid_lin.assemble_poisson(mu);
grid_lin.use_linearized_smoothers();
grid_lin.is_finest = true;

bdy     = mesh.get_boundary_node_indices(order);

n = size(K, 1);
x_gt = rand(n,1);
rhs = K*x_gt;

rhs (bdy) = 0;

% now solve and test 
maxit = min(350, size(K,1));

per = symamd(K_lin);
K_lin_chol = chol(K_lin(per,per));

% solve reordered system and revert ordering in solution
[x1,fl1,rr1,it1,rv1] = pcg(K(per,per), rhs(per), 1e-8, maxit, K_lin_chol', K_lin_chol,[]);

num_smooth = 1;
[x2,fl2,rr2,it2,rv2] = pcg(K, rhs, 1e-8, maxit, @smooth_gll);

function yo = smooth_gll(xi) 
	% first smooth 
	xs = grid_lin.smoother_chebyshev (num_smooth, xi, zeros(size(xi)));
	
  ys = xs - (K_lin \ grid_lin.residual(xi, xs));
	
	yo = grid_lin.smoother_chebyshev (num_smooth, xi, ys);
end

num_smooth = 3;
[x3,fl3,rr3,it3,rv3] = pcg(K, rhs, 1e-8, maxit, @smooth_gll);

it_gll 				= it1;
it_smooth_1 	= it2;
it_smooth_3 	= it3;

disp( [ num2str(order) ': ' num2str(it1) ' & ' num2str(it2) ' & ' num2str(it3)] );

end
