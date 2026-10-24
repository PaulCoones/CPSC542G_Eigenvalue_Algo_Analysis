%{ 
    eigenDivNConq.m provides a standalone function for computing all 
    eigenvalues of a matrix A using the divide and conquer algorithm.
    Reference: Trefethen & Bau's Numerical Linear Algebra Chapter 30
%}

function x = eigenDivNConq(A, tol)
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

function roots = secularRoots(D, beta, z, tol, k)
% Extract and sort diagonal entries
[d, idx] = sort(diag(D));
z = z(idx)

% Identify distinct poles
[poles, ~, group] = unique(d, 'stable')

% Count multiplicities
multiplicity = accumarray(group, 1)

% Combine weights associated with repeated poles
weightSquared = accumarray(group, z.^2)

% Retain groups with nonzero weights
active = weightSquared > 0

% Eigenvalues obtained directly by deflation
deflatedRoots = repelem( ...
    poles, multiplicity - double(active))

% Reduced secular problem
dReduced = poles(active);
zReduced = sqrt(weightSquared(active))

    % d = diag(D); 
    % 
    % if beta == 0
    %     warning("Vanishing rank-1 correction, " + ...
    %         "eigenvalues already in the diagonal entries of D.")
    %     roots = sort(d);
    % end
    % 
    % 
    % % Secular eqn and it's derivative
    % f = @(lambda) 1 + beta * sum(z.^2 ./ (d - lambda));
    % fp = @(lambda) beta * sum(z.^2 ./ (d - lambda).^2);
    % 
    % 
    % w = sqrt(abs(beta)) * z; 

    


end