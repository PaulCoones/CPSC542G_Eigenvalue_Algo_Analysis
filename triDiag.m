function A = triDiag(n, L)
% x = triDiag(A) Returns a symmetric tridiagonal matrix A obtained from the
% 3-point stencil central difference discretization of Poisson's equation.
% Provided an interval length L and a number of intervals n.
% Requires: n is a positive integer
%           L is positive

% Default arguments if none assigned
arguments
    n
    L = 1
end
    h = L/(n+1);
    A = (2*eye(n) - diag(ones(n-1,1),1) - diag(ones(n-1,1),-1)) / h^2;
end