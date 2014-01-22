classdef xform
  %XFORM transformation functions
  %   of the form, Xout = foo(Xin)
  
  methods(Static)
     function Xout = identity (Xin)
      Xout = Xin;
    end
 
    function Xout = twoX (Xin)
      Xout = Xin;
      Xout(:,1) = 2*Xout(:,1);
    end
    
    function Xout = twoY (Xin)
      Xout = Xin;
      Xout(:,2) = 2*Xout(:,2);
    end

		function Xout = twoSix (Xin)
			% currently only 2D 
      Xout = 2*Xin;
      Xout(:,2) = 3*Xout(:,2);
		end

    function Xout = shell (Xin)
      d = size(Xin, 2);
      R2 = 1.0; % hard coded for now.
      R1 = 0.7; % hard coded for now.
      R2byR1 = R2 / R1;
      R1sqrbyR2 = R1 * R1 / R2;
      
      if (d == 2)
        x = zeros( size(Xin(:,1)) );
        y = tan ( Xin(:,1)  * pi/4 );
        R = R1sqrbyR2 * ( R2byR1.^(Xin(:,2) + 1) ) ;
      else
        x = tan ( Xin(:,1)  * pi/4 );
        y = tan ( Xin(:,2)  * pi/4 );
        R = R1sqrbyR2 * ( R2byR1.^(Xin(:,3) + 1) );
      end
      
      q = R ./ sqrt (x.*x + y.*y + 1);
      
      if (d == 3)
        Xout(:,1) =  q.* y;
        Xout(:,2) = -q.* x;
        Xout(:,3) =  q;
      else
        Xout(:,1) =   q.* y;
        Xout(:,2) =   q;
      end
    end
  end % Static methods
  
end

% for (il = 0; il < np; ++il) {
%       /* transform abc[0] and y in-place for nicer grading */
%       x = tan (EX[il] * M_PI_4);
%       y = tan (EY[il] * M_PI_4);
% 
%       /* compute transformation ingredients */
%       R = R1sqrbyR2 * pow (R2byR1, EZ[il]);
%       q = R / sqrt (x * x + y * y + 1.);
% 
%       /* assign correct coordinates based on patch id */
%       X[il] = +q * y;
%       Y[il] = -q * x;
%       Z[il] = +q;
%     }
