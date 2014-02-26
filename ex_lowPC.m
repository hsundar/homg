
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

