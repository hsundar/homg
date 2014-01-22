classdef grid < handle
  %GRID A single grid in a multigrid heirarchy
  % compare with mgm_grid structure in mgm_multigrid.h
  
  properties
    level
    eig_max
    eig_min
    
    k_evec
    
    jacobi_omega
    jacobi_invdiag
    jacobi_inv_block_diag
    
    ssor_M
    ssor_N
    sor_omega
    
    smoother
		
    K
    L
    
    Boundary
    
    M
    
    R
    P
    
    Mesh
    
    Coarse  % handle to coarse grid 
  end % properties
  
  methods
    function grid = grid(mesh, order, coarse)
      if ((nargin < 3) || isempty(coarse))
        grid.level = 0;
        grid.Coarse = [];
      else
        grid.level = coarse.level + 1;
        grid.Coarse = coarse;
      end
      
      grid.Mesh = mesh;
      
			mesh.set_order(order); 
      grid.sor_omega = 1;
      if (~ isempty(grid.Coarse) )
         grid.P = grid.Coarse.Mesh.assemble_interpolation(order);
         grid.R = grid.P';
      end
			
			% for Dirichlet boundary conditions
			grid.Boundary = mesh.get_boundary_node_indices(order);

			%% defaults ...
		  grid.smoother = 'sor';
      grid.jacobi_omega = 2/3;
		end
    
		function assemble_poisson(grid, mu)
			% fine grid material props ...
      grid.Mesh.set_coeff (mu) ;
			% assemble for this level ... 
      [grid.K, grid.M, grid.jacobi_inv_block_diag] = grid.Mesh.assemble_poisson(grid.Mesh.order);
			syms x y z
			if ( grid.Mesh.dim == 2 )
				fx = matlabFunction(-8*pi^2*(sin(2*pi*x) * sin(2*pi*y)));
			else
				fx =matlabFunction(-12*pi^2*(sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) ));
			end
			
			grid.L = grid.Mesh.assemble_rhs(fx, grid.Mesh.order);
			grid.L(grid.Boundary) = 0;
      
      % propagate to lower grids
      if (~ isempty(grid.Coarse) )
        if isnumeric(mu)
          harmonic = 0;
          if (grid.Mesh.dim == 2)
            mu2 = reshape(mu, grid.Mesh.nelems(1), grid.Mesh.nelems(2));
            if (harmonic)
              mu_coarse = 4 ./ ( 1./mu2(1:2:end, 1:2:end) + 1./mu2(2:2:end, 1:2:end) + 1./mu2(1:2:end, 2:2:end) + 1./mu2(2:2:end, 2:2:end) );
            else
              mu_coarse = 0.25*(mu2(1:2:end, 1:2:end) + mu2(2:2:end, 1:2:end) + mu2(1:2:end, 2:2:end) + mu2(2:2:end, 2:2:end));
            end
          else
            mu3 = reshape(mu, grid.Mesh.nelems(1), grid.Mesh.nelems(2), grid.Mesh.nelems(3));
            if (harmonic)
              mu_coarse = 8 ./ ( 1./mu3(1:2:end, 1:2:end, 1:2:end) + 1./mu3(2:2:end, 1:2:end, 1:2:end) ...
                               + 1./mu3(1:2:end, 2:2:end, 1:2:end) + 1./mu3(2:2:end, 2:2:end, 1:2:end) ...
                               + 1./mu3(1:2:end, 1:2:end, 2:2:end) + 1./mu3(2:2:end, 1:2:end, 2:2:end) ...
                               + 1./mu3(1:2:end, 2:2:end, 2:2:end) + 1./mu3(2:2:end, 2:2:end, 2:2:end) );
            else
              mu_coarse = 0.125*(mu3(1:2:end, 1:2:end, 1:2:end) + mu3(2:2:end, 1:2:end, 1:2:end) + mu3(1:2:end, 2:2:end, 1:2:end) + mu3(2:2:end, 2:2:end, 1:2:end) + ...
                mu3(1:2:end, 1:2:end, 2:2:end) + mu3(2:2:end, 1:2:end, 2:2:end) + mu3(1:2:end, 2:2:end, 2:2:end) + mu3(2:2:end, 2:2:end, 2:2:end) );
            end
          end
          grid.Coarse.assemble_poisson (mu_coarse(:)) ;
        else
          grid.Coarse.assemble_poisson (mu) ;
        end
      end      
    end
    
		function set_stiffness(grid, K)
			grid.K = K;
		end
		
    % compute the residual
    function r = residual(grid, rhs, u)
      if ( nargin < 2 )
        rhs = grid.L;
      end
      if ( nargin < 3 )
        u = zeros(size(rhs));
      end

      r = grid.K*u - rhs;
    end

    function [u, rr, iter] = solve_pcg(grid, num_vcyc, smoother, v1, v2, rhs, u)
      grid.set_smoother(smoother);
      
      r = grid.residual(rhs, u);
      rho = zeros(size(u));
      
      rho = grid.vcycle(v1, v2, r, rho);
      p = rho;
      disp(['Initial residual is ' num2str(norm(r))]);
      disp('------------------------------------------');
      r0 = norm(r);
      for i=1:num_vcyc
        h = grid.K * p;
        rho_res = dot (rho, r);
        alpha = rho_res / dot ( p, h );
        u = u + alpha*p;
        r = r - alpha*h;

        disp([num2str(i, '%03d\t') ': |res| = ' num2str(norm(r),'\t%8.4e')]);
        if (norm(r)/r0 < 1e-8)
          iter = i;
          rr = norm(r)/r0;
          return;
        end
        
        % precondition ..
        rho = zeros(size(u)); 
        rho = grid.vcycle(v1, v2, r, rho);
        
        beta = dot(rho, r) / rho_res ;
        p = rho + beta*p;
      end
      disp('------------------------------------------');
      iter = num_vcyc;
      rr = norm(r)/r0;
    end
    
    function [u, rr, iter] = solve(grid, num_vcyc, smoother, v1, v2, rhs, u)
      grid.set_smoother(smoother);
      
      r = grid.residual(rhs, u);
      disp(['Initial residual is ' num2str(norm(r))]);
      disp('------------------------------------------');
      r0 = norm(r);
      for i=1:num_vcyc
        u = grid.vcycle(v1, v2, rhs, u);
        r = grid.residual(rhs, u);
        disp([num2str(i, '%03d\t') ': |res| = ' num2str(norm(r),'\t%8.4e')]);
        if (norm(r)/r0 < 1e-8)
          iter = i;
          rr = norm(r)/r0;
          return;
        end
      end
      disp('------------------------------------------');
      iter = num_vcyc;
      rr = norm(r)/r0;
    end

    % main v-cycle
    function u = vcycle(grid, v1, v2, rhs, u)
      % function u = vcycle(grid, v1, v2, rhs, u)
      % solve system using initial guess u, given rhs
      % with v1 pre and v2 post-smoothing steps
      
      % handle for the coarsest level
      if ( isempty( grid.Coarse ) )
        u = grid.K \ rhs;
        return;
      end
      
      % 1. pre-smooth
      u = grid.smooth ( v1, rhs, u );
      
      % 2. compute residual
      res = grid.residual(rhs, u);
			
      % 3. restrict
      res_coarse = grid.R * res;
      res_coarse(grid.Coarse.Boundary) = 0;
      
      % 4. recurse
      u_corr_coarse = grid.Coarse.vcycle(v1, v2, res_coarse, zeros(size(res_coarse)));
      
      % 5. prolong and correct
      u = u - grid.P * u_corr_coarse;
      
      % 6. post-smooth
      u = grid.smooth ( v2, rhs, u );
      
    end % v-cycle
    
    % smoothers
    function u = smooth (grid, v, rhs, u)
      switch(grid.smoother)
        case 'jacobi', 
          u = grid.smoother_jacobi(v, rhs, u); 
          return;
        case 'blk_jac',
          u = grid.smoother_block_jacobi(v, rhs, u); 
          return;
        case 'chebyshev', 
          u = grid.smoother_chebyshev(v, rhs, u); 
          return;
        case 'ssor',
          u = grid.smoother_sym_sor(v, rhs, u);
          return;
        otherwise
          disp('ERROR: Unrecognized smoother type'); 
          return;
      end
    end

    function set_coeff(grid, mu)
      grid.Mesh.set_coeff (mu) ;
      if (~ isempty(grid.Coarse) )
        grid.Coarse.Mesh.set_coeff (mu) ;
      end
    end
    
    function set_smoother(grid, sm)
      grid.smoother = sm;
      if (~ isempty(grid.Coarse) )
        grid.Coarse.set_smoother(sm);
      end
    end

    function u = smoother_jacobi (grid, v, rhs, u)
      % standard jacobi smoother
      if ( isempty(grid.jacobi_invdiag) )
        D = diag(grid.K);
				grid.jacobi_invdiag = 1./D;
      end
      
      for i=1:v
        res  = grid.jacobi_invdiag .* grid.residual(rhs, u);
        u = u - grid.jacobi_omega.*res;
      end  
    end % jacobi
    
    function u = smoother_block_jacobi (grid, v, rhs, u)
      % block jacobi smoother
      if ( isempty(grid.jacobi_inv_block_diag) )
        error('inv block doagonal not assembled');
      end
      
      for i=1:v
			  res  = grid.jacobi_inv_block_diag * grid.residual(rhs, u);
        u = u - grid.jacobi_omega.*res;
      end  
    end % blk jacobi
    
    function u = smoother_sym_sor (grid, v, rhs, u)
      if ( isempty ( grid.ssor_M ) )
        w = grid.sor_omega;
        n = length(u);
        grid.ssor_M = spdiags( (1/w)*diag(grid.K), 0, n, n) + tril(grid.K,-1);
        grid.ssor_N = spdiags(((1-w)/w)*diag(grid.K), 0, n, n) - triu(grid.K,1);
			end

      for i=1:v
        r = grid.residual(rhs, u);
				u = u - grid.ssor_M \ r;
        u = grid.ssor_M' \ (grid.ssor_N'*u + rhs);
      end
    end
   
    function set_sor_omega(grid, w)
      grid.sor_omega = w;
      grid.ssor_M = [];
      grid.ssor_N = [];
    end
    
		function u = smoother_chebyshev (grid, v, rhs, u)
			if ( isempty ( grid.eig_max ) )
        D = diag(grid.K);
        grid.jacobi_invdiag = 1./D;
        Kc = spdiags(grid.jacobi_invdiag,0,length(D), length(D)) * grid.K;
				opts.tol = 0.01;
				grid.eig_max = eigs(Kc, 1, 'lm', opts);  
				% grid.eig_min = eigs(Kc, 1, 'sm');  
			end

			% adjust the eigenvalues to hit the upper spectrum
			beta = grid.eig_max;
			alpha = 0.25*grid.eig_max;% (grid.eig_min + grid.eig_max)/2;

			delta = (beta - alpha)/2;
			theta = (beta + alpha)/2;
			s1 = theta/delta;
			rhok = 1./s1;

			d = zeros(size(u));

			% first loop
			res = -grid.residual ( rhs, u );
			d = res/theta.* grid.jacobi_invdiag;
			u = u + d;

			for iter = 2:v
				rhokp1 = 1/ (2*s1 - rhok);
				d1 = rhokp1 * rhok;
				d2 = 2*rhokp1 / delta;
        rhok = rhokp1;
        res = -grid.residual ( rhs, u );
				d = d1 * d + d2 * res.*grid.jacobi_invdiag;
				u = u + d;
			end
		end % chebyshev


    function [evec, eval] = get_eigenvectors(grid)
      % generate the correct matrix 
      Kc = grid.K; %(eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
      [evec, eval] = svd(full(Kc)); %eig(full(Kc), full(grid.M));
      [eval,per] = sort(diag(eval),'ascend');
      evec = evec(:,per);
    end
    

    function u0 = get_u0(grid)
      u0 = rand(size(grid.L()));
      u0(grid.Boundary) = 0;
    end
  end %methods
  
end %classdef

