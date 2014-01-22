classdef hexmesh < handle
  %HEXMESH A container class for homg regular grid meshes
  %      RGRID class for homg meshes, made of hexahedral elements.
  
  properties (SetAccess = private)
    dim=2;
    nelems=[8 8];
    order
    Xf           % transform
    
    % problem specific
    coeff
  end % properties
  
  methods
    function mesh = hexmesh(nelems, X)
      % dim = 2/3,
      % nelems = number of elements per dimension,
      % X = transform / warping function
      
      mesh.dim    = length(nelems);
      mesh.nelems  = nelems;
      mesh.Xf = X;

      mesh.coeff = @(x,y,z)(1);
    end
    
    function plot(self)
      % display the mesh. Needs X.
      figure(1);
      c = [31/256,171/256,226/256];   % default color of grid
      lw = 1;                         % default line width of grid
      
      if (self.dim == 2 )
        [x,y] = ndgrid(0:1/self.nelems(1):1.0, 0:1/self.nelems(2):1.0);
        pts = [x(:) y(:)];
        coords = self.Xf(pts);
        plot(coords(:,1), coords(:,2), 'ko');
        hold on;
        x = reshape(coords(:,1), self.nelems(1)+1, self.nelems(2)+1);
        y = reshape(coords(:,2), self.nelems(1)+1, self.nelems(2)+1);
        plot(x,y, 'Color',c,'LineWidth',lw);
        plot(x',y', 'Color',c,'LineWidth',lw);
        % boundary test
        idx = self.get_boundary_node_indices(1);
        plot(coords(idx,1), coords(idx,2), 'ro');
        title(['Quad Mesh - ' num2str(self.nelems(1)) 'x' num2str(self.nelems(2)) ])
        axis square
      else
        % 2 xy planes
        [x,y] = ndgrid( 0:1/self.nelems(1):1.0, 0:1/self.nelems(2):1.0 );
        
        % z = 0
        z = zeros (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        surf(reshape(coords(:,1), size(x)), ... 
             reshape(coords(:,2), size(x)), ...
             reshape(coords(:,3), size(x)) ...
            );
        hold on;
        x1 = reshape(coords(:,1), self.nelems(1)+1, self.nelems(2)+1);
        y1 = reshape(coords(:,2), self.nelems(1)+1, self.nelems(2)+1);
        z1 = reshape(coords(:,3), self.nelems(1)+1, self.nelems(2)+1);
        plot3(x1,y1,z1, 'k-'); 
        plot3(x1',y1',z1', 'k-'); 
        % z = 1
        z = ones  (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        surf(reshape(coords(:,1), size(x)), ... 
             reshape(coords(:,2), size(x)), ...
             reshape(coords(:,3), size(x)) ...
            );
        hold on;
        x1 = reshape(coords(:,1), self.nelems(1)+1, self.nelems(2)+1);
        y1 = reshape(coords(:,2), self.nelems(1)+1, self.nelems(2)+1);
        z1 = reshape(coords(:,3), self.nelems(1)+1, self.nelems(2)+1);
        plot3(x1,y1,z1, 'k-'); 
        plot3(x1',y1',z1', 'k-'); 
        %--------------------------------------------------------------------
        % 2 yz planes
        [y,z] = ndgrid( 0:1/self.nelems(2):1.0, 0:1/self.nelems(3):1.0 );
        % x = 0
        x = zeros (size(y));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        surf(reshape(coords(:,1), size(x)), ... 
             reshape(coords(:,2), size(x)), ...
             reshape(coords(:,3), size(x)) ...
            );
        hold on;
        x1 = reshape(coords(:,1), self.nelems(2)+1, self.nelems(3)+1);
        y1 = reshape(coords(:,2), self.nelems(2)+1, self.nelems(3)+1);
        z1 = reshape(coords(:,3), self.nelems(2)+1, self.nelems(3)+1);
        plot3(x1,y1,z1, 'k-'); 
        plot3(x1',y1',z1', 'k-'); 
        % x = 1
        x = ones (size(y));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        surf(reshape(coords(:,1), size(x)), ... 
             reshape(coords(:,2), size(x)), ...
             reshape(coords(:,3), size(x)) ...
            );
        hold on;
        x1 = reshape(coords(:,1), self.nelems(2)+1, self.nelems(3)+1);
        y1 = reshape(coords(:,2), self.nelems(2)+1, self.nelems(3)+1);
        z1 = reshape(coords(:,3), self.nelems(2)+1, self.nelems(3)+1);
        plot3(x1,y1,z1, 'k-'); 
        plot3(x1',y1',z1', 'k-'); 
        %--------------------------------------------------------------------
        % 2 xz planes
        [x,z] = ndgrid( 0:1/self.nelems(1):1.0, 0:1/self.nelems(3):1.0 );
        % y = 0
        y = zeros (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        surf(reshape(coords(:,1), size(x)), ... 
             reshape(coords(:,2), size(x)), ...
             reshape(coords(:,3), size(x)) ...
            );
        hold on;
        x1 = reshape(coords(:,1), self.nelems(1)+1, self.nelems(3)+1);
        y1 = reshape(coords(:,2), self.nelems(1)+1, self.nelems(3)+1);
        z1 = reshape(coords(:,3), self.nelems(1)+1, self.nelems(3)+1);
        plot3(x1,y1,z1, 'k-'); 
        plot3(x1',y1',z1', 'k-'); 
        % y = 1
        y = ones (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        surf(reshape(coords(:,1), size(x)), ... 
             reshape(coords(:,2), size(x)), ...
             reshape(coords(:,3), size(x)) ...
            );
        hold on;
        x1 = reshape(coords(:,1), self.nelems(1)+1, self.nelems(3)+1);
        y1 = reshape(coords(:,2), self.nelems(1)+1, self.nelems(3)+1);
        z1 = reshape(coords(:,3), self.nelems(1)+1, self.nelems(3)+1);
        plot3(x1,y1,z1, 'k-'); 
        plot3(x1',y1',z1', 'k-'); 
        
        % pretty views etc 
        view(3); axis equal;
        title(['Hex Mesh ' num2str(self.nelems(1)) 'x' num2str(self.nelems(2)) 'x' num2str(self.nelems(3))])
      end
      set (gcf, 'renderer', 'opengl');
      cameratoolbar('show');
      cameratoolbar('setmode', 'orbit');
    end
    
    function plot_fx(self, fx)
      % display the mesh. Needs X.
      hFig = figure(1);
      set(gcf,'PaperPositionMode','auto')
      set(hFig, 'Position', [200 200 800 800])
      
      if (self.dim == 2 )
        scl=8;
        [x,y] = ndgrid(0:1/(scl*self.nelems(1)):1.0, 0:1/(scl*self.nelems(2)):1.0);
        pts = [x(:) y(:)];
        coords = self.Xf(pts);
        
        ci = arrayfun( fx, coords(:,1), coords(:,2) );
        surf(reshape(coords(:,1), size(x)), ... 
             reshape(coords(:,2), size(x)), ...
             zeros(size(x)), reshape(ci, size(x)), ...
             'EdgeColor','none','LineStyle','none' ...
            );
        hold on;
        [x,y] = ndgrid(0:1/self.nelems(1):1.0, 0:1/self.nelems(2):1.0);
        pts = [x(:) y(:)];
        coords = self.Xf(pts);
        x = reshape(coords(:,1), self.nelems(1)+1, self.nelems(2)+1);
        y = reshape(coords(:,2), self.nelems(1)+1, self.nelems(2)+1);
        plot(x,y, 'k-'); 
        plot(x',y', 'k-'); 
        
        axis square;
        view(0,90);
        colorbar;
      else
        % will draw 6 planes ...
        %--------------------------------------------------------------------
        % 2 xy planes
        scl=8;
        [x,y] = ndgrid( 0:1/(scl*self.nelems(1)):1.0, 0:1/(scl*self.nelems(2)):1.0 );
        % z = 0
        z = zeros (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        ci = arrayfun( fx, coords(:,1), coords(:,2), coords(:,3) );
        surf(reshape(coords(:,1), size(x)), ... 
             reshape(coords(:,2), size(x)), ...
             reshape(coords(:,3), size(x)), ...
             reshape(ci, size(x)), ...
             'EdgeColor','none','LineStyle','none' ...
            );
        hold on;
        [x,y] = ndgrid( 0:1/self.nelems(1):1.0, 0:1/self.nelems(2):1.0 );
        z = zeros (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        x1 = reshape(coords(:,1), self.nelems(1)+1, self.nelems(2)+1);
        y1 = reshape(coords(:,2), self.nelems(1)+1, self.nelems(2)+1);
        z1 = reshape(coords(:,3), self.nelems(1)+1, self.nelems(2)+1);
        plot3(x1,y1,z1, 'k-'); 
        plot3(x1',y1',z1', 'k-'); 
        % z = 1
        [x,y] = ndgrid( 0:1/(scl*self.nelems(1)):1.0, 0:1/(scl*self.nelems(2)):1.0 );
        z = ones  (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        ci = arrayfun( fx, coords(:,1), coords(:,2), coords(:,3) );
        surf(reshape(coords(:,1), size(x)), ... 
             reshape(coords(:,2), size(x)), ...
             reshape(coords(:,3), size(x)), ...
             reshape(ci, size(x)), ...
             'EdgeColor','none','LineStyle','none' ...
            ); % 'FaceColor', 'interp', 'FaceLighting', 'phong' 
        hold on;
        [x,y] = ndgrid( 0:1/self.nelems(1):1.0, 0:1/self.nelems(2):1.0 );
        z = ones  (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        x1 = reshape(coords(:,1), self.nelems(1)+1, self.nelems(2)+1);
        y1 = reshape(coords(:,2), self.nelems(1)+1, self.nelems(2)+1);
        z1 = reshape(coords(:,3), self.nelems(1)+1, self.nelems(2)+1);
        plot3(x1,y1,z1, 'k-'); 
        plot3(x1',y1',z1', 'k-'); 
        %--------------------------------------------------------------------
        % 2 yz planes
        [y,z] = ndgrid( 0:1/(scl*self.nelems(2)):1.0, 0:1/(scl*self.nelems(3)):1.0 );
        % x = 0
        x = zeros (size(y));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        ci = arrayfun( fx, coords(:,1), coords(:,2), coords(:,3) );
        surf(reshape(coords(:,1), size(x)), ... 
             reshape(coords(:,2), size(x)), ...
             reshape(coords(:,3), size(x)), ...
             reshape(ci, size(x)), ...
             'EdgeColor','none','LineStyle','none' ...
            );
        hold on;
        [y,z] = ndgrid( 0:1/self.nelems(2):1.0, 0:1/self.nelems(3):1.0 );
        % x = 0
        x = zeros (size(y));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        x1 = reshape(coords(:,1), self.nelems(2)+1, self.nelems(3)+1);
        y1 = reshape(coords(:,2), self.nelems(2)+1, self.nelems(3)+1);
        z1 = reshape(coords(:,3), self.nelems(2)+1, self.nelems(3)+1);
        plot3(x1,y1,z1, 'k-'); 
        plot3(x1',y1',z1', 'k-'); 
        x = 1
        [y,z] = ndgrid( 0:1/(scl*self.nelems(2)):1.0, 0:1/(scl*self.nelems(3)):1.0 );
        x = ones (size(y));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        ci = arrayfun( fx, coords(:,1), coords(:,2), coords(:,3) );
        surf(reshape(coords(:,1), size(x)), ... 
             reshape(coords(:,2), size(x)), ...
             reshape(coords(:,3), size(x)), ...
             reshape(ci, size(x)), ...
             'EdgeColor','none','LineStyle','none' ...
            );
        hold on;
        [y,z] = ndgrid( 0:1/self.nelems(2):1.0, 0:1/self.nelems(3):1.0 );
        x = ones (size(y));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        x1 = reshape(coords(:,1), self.nelems(2)+1, self.nelems(3)+1);
        y1 = reshape(coords(:,2), self.nelems(2)+1, self.nelems(3)+1);
        z1 = reshape(coords(:,3), self.nelems(2)+1, self.nelems(3)+1);
        plot3(x1,y1,z1, 'k-'); 
        plot3(x1',y1',z1', 'k-'); 
        %--------------------------------------------------------------------
        % 2 xz planes
        [x,z] = ndgrid( 0:1/(scl*self.nelems(1)):1.0, 0:1/(scl*self.nelems(3)):1.0 );
        % y = 0
        y = zeros (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        ci = arrayfun( fx, coords(:,1), coords(:,2), coords(:,3) );
        surf(reshape(coords(:,1), size(x)), ... 
             reshape(coords(:,2), size(x)), ...
             reshape(coords(:,3), size(x)), ...
             reshape(ci, size(x)) ...
            );
        hold on;
        [x,z] = ndgrid( 0:1/self.nelems(1):1.0, 0:1/self.nelems(3):1.0 );
        y = zeros (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        x1 = reshape(coords(:,1), self.nelems(1)+1, self.nelems(3)+1);
        y1 = reshape(coords(:,2), self.nelems(1)+1, self.nelems(3)+1);
        z1 = reshape(coords(:,3), self.nelems(1)+1, self.nelems(3)+1);
        plot3(x1,y1,z1, 'k-'); 
        plot3(x1',y1',z1', 'k-'); 
        % y = 1
        [x,z] = ndgrid( 0:1/(scl*self.nelems(1)):1.0, 0:1/(scl*self.nelems(3)):1.0 );
        y = ones (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        ci = arrayfun( fx, coords(:,1), coords(:,2), coords(:,3) );
        surf(reshape(coords(:,1), size(x)), ... 
             reshape(coords(:,2), size(x)), ...
             reshape(coords(:,3), size(x)), ...
             reshape(ci, size(x)) ...
            );
        hold on;
        [x,z] = ndgrid( 0:1/self.nelems(1):1.0, 0:1/self.nelems(3):1.0 );
        y = ones (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        x1 = reshape(coords(:,1), self.nelems(1)+1, self.nelems(3)+1);
        y1 = reshape(coords(:,2), self.nelems(1)+1, self.nelems(3)+1);
        z1 = reshape(coords(:,3), self.nelems(1)+1, self.nelems(3)+1);
        plot3(x1,y1,z1, 'k-'); 
        plot3(x1',y1',z1', 'k-'); 
        
        % pretty views etc 
        view(3); axis square
        view(150,40);
        colorbar;
      end
      
      set (gcf, 'renderer', 'opengl');
      cameratoolbar('show');
      cameratoolbar('setmode', 'orbit');
    end
    
    function set_order(self, order)
      if isempty(self.order)
        self.order = order;
      else
        assert (order == self.order);
      end
    end
    
    function set_coeff(self, coeff)
      if ( ischar(coeff) )
        % is a string, so convert into a function
        syms x y z;
        expr = ['matlabFunction(' coeff ')'];
        self.coeff = eval(expr);
      else
        self.coeff = coeff;
      end
    end
    
    function M = assemble_mass(self, order)
      self.set_order(order);
      % assemble the mass matrix
      refel = homg.refel ( self.dim, order );
      dof = prod(self.nelems*order + 1);
      ne  = prod(self.nelems);

      % storage for indices and values
      NP = (order+1)^self.dim;
      NPNP = NP * NP;
      eM = zeros(NP, NP);
      I = zeros(ne * NPNP, 1);
      J = zeros(ne * NPNP, 1);
      val = zeros(ne * NPNP, 1);
      
      % loop over elements
      for e=1:ne
        pts = self.element_nodes(e, refel);
        detJac = self.geometric_factors(refel, pts);
        
        idx = self.get_node_indices (e, order);
        eM = self.element_mass(e, refel, detJac);
        ind1 = repmat(idx,NP,1);
        ind2 = reshape(repmat(idx',NP,1),NPNP,1);
        st = (e-1)*NPNP+1;
        en = e*NPNP;
        I(st:en) = ind1;
        J(st:en) = ind2;
        val(st:en) = eM(:);
      end
      M = sparse(I,J,val,dof,dof);
    end
    
    function K = assemble_stiffness(self, order)
      self.set_order(order);
      % assemble the stiffness matrix
      refel = homg.refel ( self.dim, order );
      
      dof = prod(self.nelems*order + 1);
      ne  = prod(self.nelems);
      
      % storage for indices and values
      NP = (order+1)^self.dim;
      NPNP = NP * NP;
      
      I = zeros(ne * NPNP, 1);
      J = zeros(ne * NPNP, 1);
      stiff_val = zeros(ne * NPNP, 1);
			
      % loop over elements
      for e=1:ne
        idx = self.get_node_indices (e, order);
        
        ind1 = repmat(idx,NP,1);
        ind2 = reshape(repmat(idx',NP,1),NPNP,1);
        st = (e-1)*NPNP+1;
        en = e*NPNP;
        I(st:en) = ind1;
        J(st:en) = ind2;
				
        pts = self.element_nodes(e, refel);
        [detJac, Jac] = self.geometric_factors(refel, pts);
        
				eMat = self.element_stiffness(e, refel, detJac, Jac);
				stiff_val(st:en) = eMat(:);
      end
    end
    
    function [K, M, iK] = assemble_poisson(self, order)
      self.set_order(order);
      % assemble the mass matrix
      refel = homg.refel ( self.dim, order );
      dof = prod(self.nelems*order + 1);
      ne  = prod(self.nelems);
      
      % storage for indices and values
      NP = (order+1)^self.dim;
      NPNP = NP * NP;
      % eMat = zeros(NP, NP);
      
      I = zeros(ne * NPNP, 1);
      J = zeros(ne * NPNP, 1);
      mass_val      = zeros(ne * NPNP, 1);
      stiff_val     = zeros(ne * NPNP, 1);
      inv_stiff_val = zeros(ne * NPNP, 1);
      ind_inner1D = repmat((2:order)', 1, order-1);
      if self.dim == 2
        ind_inner = ind_inner1D + (order+1) * (ind_inner1D'-1);
      else
        ind_inner = ind_inner1D + (order+1) * (ind_inner1D'-1);
        ind_inner = repmat(ind_inner, [1,1,order-1]);
        for i = 1:order-1
          ind_inner(:,:,i) = ind_inner(:,:,i) + i * (order+1)^2;
        end
      end

      % loop over elements
      for e=1:ne
        idx = self.get_node_indices (e, order);
        
        ind1 = repmat(idx,NP,1);
        ind2 = reshape(repmat(idx',NP,1),NPNP,1);
        st = (e-1)*NPNP+1;
        en = e*NPNP;
        I(st:en) = ind1;
        J(st:en) = ind2;
        
        pts = self.element_nodes(e, refel);
        [detJac, Jac] = self.geometric_factors(refel, pts);
        
        eMat = self.element_mass(e, refel, detJac);
        mass_val(st:en) = eMat(:);
        
        eMat = self.element_stiffness(e, refel, detJac, Jac);
        stiff_val(st:en)     = eMat(:);
      
        eMat_inner_inv = inv(eMat(ind_inner,ind_inner));
        eMat_inv = diag(diag(eMat));
        eMat_inv(ind_inner(:),ind_inner(:)) =  eMat_inner_inv;
        inv_stiff_val(st:en) = eMat_inv(:);
      end
      M = sparse(I,J,mass_val,dof,dof);
      % zero dirichlet bdy conditions
      bdy = self.get_boundary_node_indices(order);

      ii = ismember(I,bdy);
      jj = ismember(J,bdy);
      
      stiff_val = stiff_val.*(~ii).*(~jj);
      inv_stiff_val = inv_stiff_val.*(~ii).*(~jj);
      
      I = [I; bdy];
      J = [J; bdy];
      stiff_val = [stiff_val; ones(length(bdy), 1)];
      inv_stiff_val = [inv_stiff_val; ones(length(bdy), 1)];
      
      K  = sparse(I,J,stiff_val,dof,dof);
      iK = sparse(I,J,inv_stiff_val,dof,dof);
      ebdy = self.get_element_boundary_node_indices(order);
      
      iKebdry = diag(iK(ebdy,ebdy));
      iK(ebdy,ebdy) = diag(1./iKebdry);
    end

    function f = assemble_rhs(self, fx, order)
      self.set_order(order);
      refel = homg.refel ( self.dim, order );
      
      dof = prod(self.nelems*order + 1);
      ne  = prod(self.nelems);
      
      f = zeros(dof,1);
      
      % loop over elements
      for e=1:ne
        idx = self.get_node_indices (e, order);
        
        pts = self.element_nodes(e, refel);
        J = self.geometric_factors(refel, pts);
        
        gpts =  self.element_gauss(e, refel);
        
        if (self.dim == 2)
          fd = arrayfun( fx, gpts(:,1), gpts(:,2) );
        else
          fd = arrayfun( fx, gpts(:,1), gpts(:,2), gpts(:,3) );
        end
        
        Jd = refel.W .* J .* fd; 
        
        f(idx) = f(idx) + refel.Q' * Jd;
      end
    end
    
    function P = assemble_interpolation(self, order)
      % assemble prolongation operator from coarse (self) to fine mesh
			refel = homg.refel ( self.dim, self.order );
      
      if ( order == self.order )
				dof_coarse = prod(  self.nelems * self.order + 1);
        dof_fine   = prod(2*self.nelems * self.order + 1);
        NP_c = (self.order+1)^self.dim;
				NP_f = (2*self.order+1)^self.dim;
        Pe = refel.Ph;
      else
        assert (order == 2*self.order);
        NP_c = (self.order+1)^self.dim;
				NP_f = (order+1)^self.dim;
        dof_coarse = prod(self.nelems * self.order + 1);
        dof_fine   = prod(self.nelems * order + 1);
				Pe = refel.Pp; 
			end
			
      ne  = prod(self.nelems);

      % storage for indices and values
      NPNP = NP_c * NP_f;

      I = zeros(ne * NPNP, 1);
      J = zeros(ne * NPNP, 1);
      val = zeros(ne * NPNP, 1);

      for e=1:ne
        [idx_c, idx_f] = self.get_interpolation_indices (e);

        ind1 = repmat(idx_f,NP_c,1);
        ind2 = reshape(repmat(idx_c',NP_f,1),NPNP,1);
        st = (e-1)*NPNP+1;
        en = e*NPNP;
        I(st:en) = ind1;
        J(st:en) = ind2;

        val(st:en) = Pe(:);
      end

      [u_ij,q] = unique([I,J],'rows','first');
      u_val   = val(q);
      I = u_ij(:,1);
      J = u_ij(:,2);

      P = sparse (I,J,u_val,dof_fine,dof_coarse);
    end

    function idx = get_node_indices ( self, eid, order )
    % determine global node indices for a given element
    if ( self.dim == 2)
        [i,j] = ind2sub (self.nelems, eid);

        i_low   = (i-1)*order + 1;   i_high =  i*order + 1;
        j_low   = (j-1)*order + 1;   j_high =  j*order + 1;

        [i,j] = ndgrid(i_low:i_high, j_low:j_high);

        idx     = sub2ind (self.nelems*order + 1, i(:), j(:));
      else
        [i,j,k] = ind2sub (self.nelems, eid);

        i_low   = (i-1)*order + 1;   i_high =  i*order + 1;
        j_low   = (j-1)*order + 1;   j_high =  j*order + 1;
        k_low   = (k-1)*order + 1;   k_high =  k*order + 1;

        [i,j,k] = ndgrid(i_low:i_high, j_low:j_high, k_low:k_high);

        idx     = sub2ind (self.nelems*order + 1, i(:), j(:), k(:) );
      end
    end

    function idx = get_linear_node_indices ( self, eid, order )
    % determine global node indices for a given element
    if ( self.dim == 2)
        [i,j] = ind2sub (self.nelems*order, eid);

        [i,j] = ndgrid(i:i+1, j:j+1);

        idx     = sub2ind (self.nelems*order + 1, i(:), j(:));
      else
        [i,j,k] = ind2sub (self.nelems*order, eid);

        [i,j,k] = ndgrid(i:i+1, j:j+1, k:k+1);

        idx     = sub2ind (self.nelems*order + 1, i(:), j(:), k(:) );
      end
    end

    function [idx_coarse, idx_fine] = get_interpolation_indices ( self, eid )
    % determine global node indices for a given element
    if ( self.dim == 2)
        [i,j] = ind2sub (self.nelems, eid);

        i_low       = (i-1)*self.order + 1;   i_high =  i*self.order + 1;
        j_low       = (j-1)*self.order + 1;   j_high =  j*self.order + 1;
        [i,j]       = ndgrid(i_low:i_high, j_low:j_high);
        idx_coarse  = sub2ind (self.nelems*self.order + 1, i(:), j(:));

        [i,j]       = ndgrid(2*i_low-1:2*i_high-1, 2*j_low-1:2*j_high-1);
        idx_fine    = sub2ind (2*self.nelems*self.order + 1, i(:), j(:));
      else
        [i,j,k] = ind2sub (self.nelems, eid);

        i_low       = (i-1)*self.order + 1;   i_high =  i*self.order + 1;
        j_low       = (j-1)*self.order + 1;   j_high =  j*self.order + 1;
        k_low       = (k-1)*self.order + 1;   k_high =  k*self.order + 1;
        [i,j,k]     = ndgrid(i_low:i_high, j_low:j_high, k_low:k_high);
        idx_coarse  = sub2ind (self.nelems*self.order + 1, i(:), j(:), k(:) );

        [i,j,k]     = ndgrid(2*i_low-1:2*i_high-1, 2*j_low-1:2*j_high-1, 2*k_low-1:2*k_high-1);
        idx_fine    = sub2ind (2*self.nelems*self.order + 1, i(:), j(:), k(:) );
      end
    end

    function idx = get_boundary_node_indices(self, order)
    % function idx = get_boundary_node_indices(self, order)
      %    returns indices of boundary nodes, for setting 
      %    boundary conditions       
      if (self.dim == 2)
        [x,y] = ndgrid(1:self.nelems(1)*order+1,1:self.nelems(2)*order+1);
        
        idx = [ find(x == 1);     find(x == (self.nelems(1)*order+1));
                find(y == 1);     find(y == (self.nelems(2)*order+1))  ];
        
        idx = unique(idx);
      else 
         [x,y,z] = ndgrid(1:self.nelems(1)*order+1,1:self.nelems(2)*order+1,1:self.nelems(3)*order+1);
         
         idx = [ find(x == 1);     find(x == (self.nelems(1)*order+1));
                 find(y == 1);     find(y == (self.nelems(2)*order+1));
                 find(z == 1);     find(z == (self.nelems(3)*order+1))  ];
        
        idx = unique(idx);
      end
        end
        
		function idx = get_element_boundary_node_indices(self, order)
	  % function idx = get_element_boundary_node_indices(self, order)
      %    returns indices of element boundary nodes, for block 
      %    Jacobi smoother       
      if (self.dim == 2)
        [x,y] = ndgrid(1:self.nelems(1)*order+1,1:self.nelems(2)*order+1);
        
        idx = [ find(mod(x,order) == 1);
                find(mod(y,order) == 1);];
        
        idx = unique(idx);
      else 
         [x,y,z] = ndgrid(1:self.nelems(1)*order+1,1:self.nelems(2)*order+1,1:self.nelems(3)*order+1);
         
         idx = [ find(mod(x,order) == 1);
                 find(mod(y,order) == 1);
                 find(mod(z,order) == 1);];
        
        idx = unique(idx);
      end
        end
    
    function Me = element_mass(self, eid, refel, J)
      % element mass matrix
      Md = refel.W .* J ; 
      
      Me = refel.Q' * diag(Md) * refel.Q;
    end
    
    function Ke = element_stiffness(self, eid, r, J, D)
      % element mass matrix
      
%             | Qx Qy Qz || rx ry rz |     | rx sx tx || Qx |
%    Ke =                 | sx sy sz | J W | ry sy ty || Qy |
%                         | tx ty tz |     | rz sz tz || Qz |

     	gpts = self.element_gauss(eid, r);

			nn = length(J);
			
      factor = zeros(nn, 6);

      %             1  4  5
      % factor      4  2  6
      %             5  6  3
      
      if (self.dim == 2 )
        mu = arrayfun( self.coeff, gpts(:,1), gpts(:,2) );
        
        factor (:,1) = (D.rx.*D.rx + D.ry.*D.ry ) .* J .* r.W .* mu ; % d2u/dx^2
        factor (:,2) = (D.sx.*D.sx + D.sy.*D.sy ) .* J .* r.W .* mu ; % d2u/dy^2
        factor (:,3) = (D.rx.*D.sx + D.ry.*D.sy ) .* J .* r.W .* mu ; % d2u/dxdy
        
        Ke =   r.Qx' * diag(factor(:,1)) * r.Qx ...
             + r.Qy' * diag(factor(:,2)) * r.Qy ...
             + r.Qx' * diag(factor(:,3)) * r.Qy ...
             + r.Qy' * diag(factor(:,3)) * r.Qx ;
      else
        mu = arrayfun( self.coeff, gpts(:,1), gpts(:,2), gpts(:,3) );
				
        % first compute dj.w.J.J'
        factor (:,1) = (D.rx.*D.rx + D.ry.*D.ry + D.rz.*D.rz ) .* J .* r.W .* mu ; % d2u/dx^2
        factor (:,2) = (D.sx.*D.sx + D.sy.*D.sy + D.sz.*D.sz ) .* J .* r.W .* mu ; % d2u/dy^2
        factor (:,3) = (D.tx.*D.tx + D.ty.*D.ty + D.tz.*D.tz ) .* J .* r.W .* mu ; % d2u/dz^2
        
        factor (:,4) = (D.rx.*D.sx + D.ry.*D.sy + D.rz.*D.sz ) .* J .* r.W .* mu ; % d2u/dxdy
        factor (:,5) = (D.rx.*D.tx + D.ry.*D.ty + D.rz.*D.tz ) .* J .* r.W .* mu ; % d2u/dxdz
        factor (:,6) = (D.sx.*D.tx + D.sy.*D.ty + D.sz.*D.tz ) .* J .* r.W .* mu ; % d2u/dydz
        
        Ke =   r.Qx' * diag(factor(:,1)) * r.Qx ...
             + r.Qy' * diag(factor(:,2)) * r.Qy ...
             + r.Qz' * diag(factor(:,3)) * r.Qz ...
             + r.Qx' * diag(factor(:,4)) * r.Qy ...
             + r.Qy' * diag(factor(:,4)) * r.Qx ...
             + r.Qx' * diag(factor(:,5)) * r.Qz ...
             + r.Qz' * diag(factor(:,5)) * r.Qx ...
             + r.Qz' * diag(factor(:,6)) * r.Qy ...
             + r.Qy' * diag(factor(:,6)) * r.Qz ;
      end
      
    end
    
    function [J, D] = geometric_factors( self, refel, pts )
      % change to using Qx etc ?
      if (refel.dim == 1)
        xr  = refel.Dg * pts; 
        J = xr; 
      elseif (refel.dim == 2)
        [xr, xs] = homg.tensor.grad2 (refel.Dg, pts(:,1));
        [yr, ys] = homg.tensor.grad2 (refel.Dg, pts(:,2));
        
        J = -xs.*yr + xr.*ys;
      else
        [xr, xs, xt] = homg.tensor.grad3 (refel.Dg, pts(:,1));
        [yr, ys, yt] = homg.tensor.grad3 (refel.Dg, pts(:,2));
        [zr, zs, zt] = homg.tensor.grad3 (refel.Dg, pts(:,3));
        
        J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
      end
      
      if (nargout > 1)
        if (refel.dim == 1)
          D.rx = 1./J;
        elseif (refel.dim == 2)
          D.rx =  ys./J;
          D.sx = -yr./J;
          D.ry = -xs./J;
          D.sy =  xr./J;
        else
          D.rx =  (ys.*zt - zs.*yt)./J;
          D.ry = -(xs.*zt - zs.*xt)./J;
          D.rz =  (xs.*yt - ys.*xt)./J;
          
          D.sx = -(yr.*zt - zr.*yt)./J;
          D.sy =  (xr.*zt - zr.*xt)./J;
          D.sz = -(xr.*yt - yr.*xt)./J;
          
          D.tx =  (yr.*zs - zr.*ys)./J;
          D.ty = -(xr.*zs - zr.*xs)./J;
          D.tz =  (xr.*ys - yr.*xs)./J;
        end
        
      end
    end
    
    function coords = linear_element_nodes(self, elem, order) 
      
      if (self.dim == 2)
        [i,j] = ind2sub (self.nelems*order, elem);
        
        x1d = homg.hexmesh.getGLLcoords(order, self.nelems(1));
        y1d = homg.hexmesh.getGLLcoords(order, self.nelems(2));
        
        [x, y] = ndgrid(x1d(i:i+1), y1d(j:j+1));
        pts = [x(:) y(:)];
      else
        [i,j,k] = ind2sub (self.nelems*order, elem);
        
        x1d = homg.hexmesh.getGLLcoords(order, self.nelems(1));
        y1d = homg.hexmesh.getGLLcoords(order, self.nelems(2));
        z1d = homg.hexmesh.getGLLcoords(order, self.nelems(3));
        
        [x, y, z] = ndgrid(x1d(i:i+1), y1d(j:j+1), z1d(k:k+1));
        pts = [x(:) y(:) z(:)];
      end
      
      coords = self.Xf(pts);
    end
    
    function coords = element_nodes(self, elem, refel)
      h = 1./self.nelems;
      
      if ( self.dim == 2)
        [i,j] = ind2sub (self.nelems, elem);
        idx = [i j];
      else
        [i,j,k] = ind2sub (self.nelems, elem);
        idx = [i j k];
      end
      
      p_mid = (idx - 0.5) .* h;
      p_gll = refel.r * 0.5 * h;
      nodes = bsxfun(@plus, p_mid, p_gll) ;
      
      if ( self.dim == 2)
        [x, y] = ndgrid(nodes(:,1), nodes(:,2));
        pts = [x(:) y(:)];
      else
        [x, y, z] = ndgrid(nodes(:,1), nodes(:,2), nodes(:,3));
        pts = [x(:) y(:) z(:)];
      end
      
      coords = self.Xf(pts);
    end
		
    function coords = element_gauss(self, elem, refel)
      % function pts = element_gauss(self, elem, refel)
      % returns location of gauss coordinates of order
      % for element 
      
      if (self.order == refel.N) 
        h = 1./self.nelems;
      
        if ( self.dim == 2)
          [i,j] = ind2sub (self.nelems, elem);
          idx = [i j];
        else
          [i,j,k] = ind2sub (self.nelems, elem);
          idx = [i j k];
        end
      
        p_mid = (idx - 0.5) .* h;
        p_gau = refel.g * 0.5 * h;
        nodes = bsxfun(@plus, p_mid, p_gau) ;
        
        if ( self.dim == 2)
          [x, y] = ndgrid(nodes(:,1), nodes(:,2));
          pts = [x(:) y(:)];
        else
          [x, y, z] = ndgrid(nodes(:,1), nodes(:,2), nodes(:,3));
          pts = [x(:) y(:) z(:)];
        end
      else
        assert(refel.N == 1); 
        % ... get gll points ...
        if (self.dim == 2)
          [i,j] = ind2sub (self.nelems*self.order, elem);

          x1d = homg.hexmesh.getGLLcoords(self.order, self.nelems(1));
          y1d = homg.hexmesh.getGLLcoords(self.order, self.nelems(2));

          xg = x1d(i) + (x1d(i+1) - x1d(i))*(refel.g + 1)*0.5;
          yg = y1d(j) + (y1d(j+1) - y1d(j))*(refel.g + 1)*0.5;
          
          [x, y] = ndgrid(xg, yg);
          
          pts = [x(:) y(:)];
        else
          [i,j,k] = ind2sub (self.nelems*self.order, elem);

          x1d = homg.hexmesh.getGLLcoords(self.order, self.nelems(1));
          y1d = homg.hexmesh.getGLLcoords(self.order, self.nelems(2));
          z1d = homg.hexmesh.getGLLcoords(self.order, self.nelems(3));
  
          xg = x1d(i) + (x1d(i+1) - x1d(i))*(refel.g + 1)*0.5;
          yg = y1d(j) + (y1d(j+1) - y1d(j))*(refel.g + 1)*0.5;
          zg = z1d(k) + (z1d(k+1) - z1d(k))*(refel.g + 1)*0.5;
          
          [x, y, z] = ndgrid(xg, yg, zg);
          
          pts = [x(:) y(:) z(:)];
        end
      end
      
      coords = self.Xf(pts);
    end
    
    function [K, M] = assemble_poisson_linearized (self, order)
      self.set_order(order);
      
      refel = homg.refel ( self.dim, 1 );
      
      dof = prod ( self.nelems*order + 1);
      ne  = prod ( self.nelems*order ) ;
      
      % storage for indices and values
      NP = (1+1)^self.dim; % linear elements
      NPNP = NP * NP;
      % eMat = zeros(NP, NP);
      
      I = zeros(ne * NPNP, 1);
      J = zeros(ne * NPNP, 1);
      mass_val = zeros(ne * NPNP, 1);
      stiff_val = zeros(ne * NPNP, 1);
      
      % loop over elements
      for e=1:ne
        idx = self.get_linear_node_indices (e, order);
        
        ind1 = repmat(idx,NP,1);
        ind2 = reshape(repmat(idx',NP,1),NPNP,1);
        st = (e-1)*NPNP+1;
        en = e*NPNP;
        I(st:en) = ind1;
        J(st:en) = ind2;
        
        pts = self.linear_element_nodes(e, order);
        
        [detJac, Jac] = self.geometric_factors(refel, pts);
        
        eMat = self.element_mass(e, refel, detJac);
        mass_val(st:en) = eMat(:);
        
        eMat = self.element_stiffness(e, refel, detJac, Jac);
        stiff_val(st:en) = eMat(:);
      end
      M = sparse(I,J,mass_val,dof,dof);
      % zero dirichlet bdy conditions
      bdy = self.get_boundary_node_indices(order);

      ii = ismember(I,bdy);
      jj = ismember(J,bdy);
      
      stiff_val = stiff_val.*(~ii).*(~jj);
      I = [I; bdy];
      J = [J; bdy];
      stiff_val = [stiff_val; ones(length(bdy), 1)];
      
      K = sparse(I,J,stiff_val,dof,dof);
    end
 
  end % methods
  
  methods(Static) 
    function coords = getGLLcoords(order, elems)
      % function coords=getGLLcoords(order, elems)
      % returns location of gll coordinates of order
      % for elements in [0,1]
      
      fac = 1.0/(2*elems);
      
      % gll coordinates in [-1,1]
      x = homg.basis.gll (0,0,order)';
      
      x = (x + 1)*fac;
      
      
      coords = [];
      for i=1:elems
        y = x + (i-1)/elems;
        coords = [coords y(1:end-1)];
      end
      
      coords = [coords 1.0];
    end
      
    function coords = getUniformCoords(order, elems)
      coords = linspace(0, 1, order*elems+1);
    end
		
		function coords = getElementCenters(order, elems)
			% order is ignored ... 
			nodes = linspace(0,1, elems+1);
			coords = 1/2*(nodes(1:end-1) + nodes(2:end));
		end
      
    function C = stats(nelems, order)
      % function Ch = stats(nelems, order)
      %   given number of elements and the order,
      %   this function calculates different node 
      %   stats for the mesh
      d               = length(nelems);
      C.num_nodes     = prod(nelems*order + 1);
      C.num_elements  = prod(nelems);
      C.num_bdy_nodes = C.num_nodes - prod(nelems*order - 1);
      
      C.num_int_elements = prod(nelems - 1);
      C.num_bdy_elements = C.num_elements - C.num_int_elements;
      
      C.nnz = (order+2)^d*C.num_nodes;
%       if (d == 2)
%         
%       else
%         
%       end
    end
  end  % static methods
  
  
end % class
