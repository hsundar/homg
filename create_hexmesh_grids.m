function grid = create_hexmesh_grids(dim, mu, xform, orders, nelems)
% function grid = create_hexmesh_grids(dim, mu, xform, orders, nelems)
%
%   dim  is the dimension of the mesh, either 2 or 3
%   mu is a function for the coefficients, of the form @(x,y) or @(x,y,z)
%   xform is the transformation or warping, examples is +homg.xform
%   
%   specify nelems as an array of sizes, preferably factors of 2
%     from coarse to fine. smallest size first
%   specify orders as an array of sizes, as multiples of 2, and 
%     going down to 1. lowest order first
%
%   returns handle to the finest grid.
%
% example:
%          mu = @(x,y)(1 + 1e6*( cos(2*pi*x)^2 + cos(2*pi*y)^2 ) );
%          fine_grid = create_grid_heirarchy(2, mu, @homg.xform.identity, [1 2 4], [4 8 16]);
%
% @author: Hari Sundar - hari@ices.utexas.edu 
% @date  : 28-Aug-2012, 31-Apr-2013

num_hgrids = length(nelems);
num_pgrids = length(orders);

num_grids = num_hgrids + num_pgrids - 1;

m = homg.hexmesh(repmat(nelems(1), 1, dim), xform);
coarse = homg.grid(m, orders(1));

% create h-grids
for i=2:num_hgrids
  m = homg.hexmesh(repmat(nelems(i), 1, dim), xform);
  grid = homg.grid(m, orders(1), coarse);
  coarse = grid;
end

hfine = nelems(num_hgrids);

% create p-grids
for i=2:num_pgrids
  m = homg.hexmesh(repmat(hfine, 1, dim), xform);
  grid = homg.grid(m, orders(i), coarse);
  coarse = grid;
end

% Now assemble matrices ...
grid.assemble_poisson(mu);

%% test for exact galerkin operator ...  
% coarse = grid;
% while ( ~ isempty(coarse.Coarse) )
% 	RKP = coarse.R * coarse.K * coarse.P;
% 	coarse.Coarse.set_stiffness(RKP);
% 	coarse = coarse.Coarse;	
% end
