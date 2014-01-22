classdef tensor
    %TENSOR Basic tensor routines
    
    methods(Static)
        % 2D routines
        function y = IAX (A, x)
            N = size (A, 1);
            y = A * reshape(x, N, N);
            y = y(:);
        end
        
        function y = AIX (A, x)
            N = size (A, 1);
            y = A * reshape(x, N, N)';
            y = y'; 
            y = y(:);
        end
        
        % 3D routines
        function y = IIAX (A, x)
            N = size (A, 1);
            y = A * reshape(x, N, N*N);
            y = y(:);
        end
        
        function y = IAIX (A, x)
            N = size (A, 1);
            q = reshape(x, N, N, N);
            y = zeros(N,N,N);
            for i=1:N
                y(i,:,:) = A * squeeze( q(i,:,:) );
            end
            y = y(:);
        end
        
        function y = AIIX (A, x)
            N = size (A, 1);
            y = reshape(x, N*N, N) * A';
            y = y(:);
        end
        
        function du = grad(refel, u)
            du = zeros(length(u), refel.dim);
            if (refel.dim == 2)
                du(:,1) = homg.tensor.IAX (refel.Dr, u);
                du(:,2) = homg.tensor.AIX (refel.Dr, u);
            else
                du(:,1) = homg.tensor.IIAX (refel.Dr, u);
                du(:,2) = homg.tensor.IAIX (refel.Dr, u);
                du(:,3) = homg.tensor.AIIX (refel.Dr, u);
            end
        end
        
        function [dx, dy] = grad2(A, x)
           dx = homg.tensor.IAX (A, x);
           dy = homg.tensor.AIX (A, x);
        end
        
        function [dx, dy, dz] = grad3(A, x)
           dx = homg.tensor.IIAX (A, x);
           dy = homg.tensor.IAIX (A, x);
           dz = homg.tensor.AIIX (A, x);
        end
    end
    
end

