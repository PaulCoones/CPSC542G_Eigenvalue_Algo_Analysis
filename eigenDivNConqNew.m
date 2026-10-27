%{ 
    eigenDivNConq.m provides a standalone function for computing all 
    eigenvalues of a matrix A using the divide and conquer algorithm.
    Reference: Trefethen & Bau's Numerical Linear Algebra Chapter 30
%}

function [lambda, first, last] = eigenDivNConqNew(A, tol, maxIters)
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

arguments
    A
    tol = 1e-8
    maxIters = 200
end

m = size(A, 1);

% Base case : 1x1 matrix, eigenvalues and vectors are trivial
% Future consideration: Try a 2x2 matrix base case!
if m == 1
    lambda = A(1,1);
    first = 1;
    last = 1;
    return;
end

% Divide the subproblem
n = floor(m/2);
beta = A(n, n+1);
T1 = A(1:n, 1:n);
T2 = A(n+1:m, n+1:m);

% Adjust for the rank-1 update matrix
T1(end, end) = T1(end, end) - beta;
T2(1,1) = T2(1,1) - beta;


% Work on the divided ~halves
[lambda1, first1, last1] = eigenDivNConqNew(T1, tol);
[lambda2, first2, last2] = eigenDivNConqNew(T2, tol);

% Construct diagonal matrix and rank-one vector
D = diag([lambda1; lambda2]);
z = [last1'; first2'];

% Construct boundary rows inherited from the recursive subproblems
firstRowBound = [first1, zeros(1,m-n)];
lastRowBound = [zeros(1,n), last2];

% Return eigenvalues and boundary rows
[lambda, first, last] = secularBoundary(D, beta, z, firstRowBound, ...
    lastRowBound, tol, maxIters);

return;