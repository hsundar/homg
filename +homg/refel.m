classdef refel < handle
    %REFEL hexahedral reference element
    %    - Equivalent to mangll_refel_t
    
    properties
        dim
        N      % polynomial order
        Nrp    % number of 1D interpolation points
        
        r      % 1D reference coordinates of the interpolation nodes ( 1 x Nrp )
        g      % 1D reference coordinates of the gauss quadrature points
        w      % 1d weights for gauss quadrature
        wgll   % 1d weights for gll   quadrature
        
        W      % dim-dimensional weights (gauss) 
        Wgll   % dim-dimensional weights (gll)
        
        Wfgll  % face weights (gll)
        
        Vr     % 1D Vandermonde matrix of Legendre polynomials at r 
        gradVr %     and their derivative (Nrp x Nrp)
        
        Vg     % 1D Vandermonde matrix of Legendre polynomials at g 
        gradVg %     and their derivative (Nrp x Nrp)
        
        Dr     % Derivative of Lagrange interpolants at the interpolation nodes.
               % ( Nrp x Nrp )
               % Dr(i,j) = lagrange_i' (r_j)
        Dg     % Derivative of Lagrange interpolants at the gauss quad points
               % ( Nrp x Nrp )
               % Dr(i,j) = lagrange_i' (r_j)
        
        Q      % map to gauss points 
    
        Qx
        Qy
        Qz
        
        
        % Prolongation 
        Ph      % interpolation from this element to its 4/8 children
        Pp      % interpolation from this element to its 2p version
				
        Mr     % exact 1D Mass matrix (Nrp x Nrp) at gll
        invMr  % and its inverse
        
        
				Mg     % exact 1D Mass matrix (Nrp x Nrp) at gauss
				invMg
    end
    
    methods
        function elem = refel(d, order)
            % Setup a d-dimensional reference element 
            % order = polynomial order of the reference element
            elem.dim    = d;
            elem.N      = order;
            elem.Nrp    = order + 1;
            
            [elem.r, elem.wgll]  = homg.basis.gll (0, 0, elem.N);
            
            r_hby2      = [0.5*(elem.r - 1); 0.5*(elem.r(2:end) + 1)];
            r_2p        = homg.basis.gll (0, 0, 2*elem.N);
						
            [elem.g, elem.w] = homg.basis.gauss(0, 0, elem.N);
            
            elem.Vr     = zeros (order+1, order+1);
            elem.gradVr = zeros (order+1, order+1);
            
            elem.Vg     = zeros (order+1, order+1);
            elem.gradVg = zeros (order+1, order+1);
            
            Vph     = zeros (order+1, 2*order+1);
						Vpp     = zeros (order+1, 2*order+1);
            
            for i=1:elem.Nrp
                elem.Vr(i,:)     = homg.basis.polynomial (elem.r, 0, 0, i-1);
                elem.gradVr(i,:) = homg.basis.gradient (elem.r, 0, 0, i-1);
                
                elem.Vg(i,:)     = homg.basis.polynomial (elem.g, 0, 0, i-1);
                elem.gradVg(i,:) = homg.basis.gradient (elem.g, 0, 0, i-1);
                
                Vph(i,:)          = homg.basis.polynomial (r_hby2, 0, 0, i-1);
								Vpp(i,:)          = homg.basis.polynomial (r_2p,   0, 0, i-1);
            end
        
            elem.Dr     = transpose(elem.Vr \ elem.gradVr);
            
            elem.Dg     = transpose(elem.Vr \ elem.gradVg);
            
            iVr         = elem.Vr \ eye(order+1);
						iVg         = elem.Vg \ eye(order+1);
            
            q1d         = transpose (elem.Vr \ elem.Vg);  
            
						p_h_1d      = transpose (elem.Vr \ Vph);  
            p_p_1d      = transpose (elem.Vr \ Vpp);  
						
            elem.W      = zeros(elem.Nrp^elem.dim, 1);
            elem.Wgll   = zeros(elem.Nrp^elem.dim, 1);
            if (d == 2)
              elem.Q  = kron(q1d, q1d) ;
              
							elem.Ph = kron(p_h_1d, p_h_1d) ;
							elem.Pp = kron(p_p_1d, p_p_1d) ;
              
              elem.Qx = kron(q1d, elem.Dg);
              elem.Qy = kron(elem.Dg, q1d);
              
              sk = 1;
              for i=1:elem.Nrp
                for j=1:elem.Nrp
                  elem.W(sk)    = elem.w(i) * elem.w(j);
                  elem.Wgll(sk) = elem.wgll(i) * elem.wgll(j);
                  sk = sk + 1;
                end
              end
              
              Wfgll = elem.wgll;
            else
              elem.Q  = kron(kron(q1d, q1d), q1d);
							
              elem.Ph = kron(kron(p_h_1d, p_h_1d), p_h_1d);
							elem.Pp = kron(kron(p_p_1d, p_p_1d), p_p_1d);
              
              elem.Qx = kron(kron(q1d, q1d), elem.Dg);
              elem.Qy = kron(kron(q1d, elem.Dg), q1d);
              elem.Qz = kron(kron(elem.Dg, q1d), q1d);
              
              sk = 1;
              sf = 1;
              for i=1:elem.Nrp
                for j=1:elem.Nrp
                  for k=1:elem.Nrp
                    elem.W(sk)    = elem.w(i) * elem.w(j) * elem.w(k);
                    elem.Wgll(sk) = elem.wgll(i) * elem.wgll(j) * elem.wgll(k);
                    sk = sk + 1;
                  end
                  elem.Wfgll(sf) = elem.wgll(i) * elem.wgll(j);
                  sf = sf + 1;
                end
              end
            end
            
            elem.Mr     = iVr * iVr';
            elem.invMr  = elem.Mr \ eye(order+1);
            elem.Mg     = iVg * iVg';
            elem.invMg  = elem.Mg \ eye(order+1);
            
        end
    end
    
end

