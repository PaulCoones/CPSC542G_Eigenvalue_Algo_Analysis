%{ 
    eigenDivNConq.m provides a standalone function for computing all 
    eigenvalues of a matrix A using the divide and conquer algorithm.
    Reference: Trefethen & Bau's Numerical Linear Algebra Chapter 30
%}

function [lambda, first, last] = eigenDivNConq(A, tol)
% x = eigenDivNConq(A) Returns the eigenvalues of A in increasing order, 
%   by applying the divide and conquer eigenvalue algorithm, that is by 
%   recursively subdividing the problem into problems of smaller dimension.
% Requires: A is a symmetric tridiagonal square nxn matrix with n >= 1
%           A is irreducible s.t. all its off-diagonal entries are nonzero
%           tol is a finite positive number and ideally made small
% Time Complexity: ?? ADD ME ??
%                  Where n is number of rows/cols in the matrix A
%                  U and L are the upper and lower bounds on the
%                  eigenvalues from Gershgorin's circle theorem
%                  tol is the interval length tolerance
% Example Input: ?? ADD ME ??
% Example Output: ?? ADD ME ??

m = size(A, 1);

% Ensure we have some work to do
if m == 1
    lambda = A(1,1);
    first = 1;
    last = 1;
    return;
end

n = floor(m/2);
beta = A(n, n+1);
T1 = A(1:n, 1:n);
T2 = A(n+1:m, n+1:m);
T1(end, end) = T1(end, end) - beta;
T2(1,1) = T2(1,1) - beta;

[lambda1, first1, last1] = eigenDivNConq(T1, tol);
[lambda2, first2, last2] = eigenDivNConq(T2, tol);

% Construct diagonal matrix and rank-one vector
D = diag([lambda1; lambda2])
z = [last1'; first2']

% Solve secular equation
[lambda, V] = secularRoots(D, beta, z, tol, 200)

% Propagate only the boundary eigenvector components
first = [first1, zeros(1,m-n)] * V
last = [zeros(1,n), last2] * V

return;

% m = floor(n/2);
% beta = A(m,m+1);
% T1 = A(1:m,1:m);
% T2 = A(m+1:n,m+1:n);
% T1(end,end) = T1(end,end) - beta;
% T2(1,1) = T2(1,1) - beta;
% 

% Algorithm: divideConquer(A)
% 
% n = size(A,1)
% 
% IF n == 1
% return A(1,1), 1, 1
% END
% 
% m = floor(n/2)
% beta = A(m,m+1)
% 
% A1 = A(1:m,1:m)
% A2 = A(m+1:n,m+1:n)
% 
% A1(m,m) = A1(m,m) - beta
% A2(1,1) = A2(1,1) - beta
% 
% % Recursively solve subproblems
% [lambda1, first1, last1] = divideConquer(A1)
% [lambda2, first2, last2] = divideConquer(A2)
% 
% % Construct diagonal matrix and rank-one vector
% D = diag([lambda1; lambda2])
% z = [last1'; first2']
% 
% % Solve secular equation
% [lambda, V] = secularEigenpairs(D, beta, z)
% 
% % Propagate only the boundary eigenvector components
% first = [first1, zeros(1,n-m)] * V
% last = [zeros(1,m), last2] * V
% 
% return lambda, first, last

end