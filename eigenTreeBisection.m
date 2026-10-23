%{ 
    eigenTreeBisection.m provides a standalone function for computing all 
    eigenvalues of a matrix A using the Bisection eigenvalue algorithm 
    sequentially and intelligently using an interval decision tree to avoid
    repeated calculations like in eigenBisection.m's naive implementation.
    Reference: Trefethen & Bau's Numerical Linear Algebra Chapter 30
%}

function x = eigenTreeBisection(A, tol)
% x = eigenTreeBisection(A) Returns the eigenvalues of A in increasing order, 
%   by applying the Bisection eigenvalue algorithm for each eigenvalue 
%   until an interval length tolerance tol is reached.
%   This implementation is sequential, and is not trivially parallelizable
%   due to the way in which we split intervals.
% Requires: A is a symmetric tridiagonal square nxn matrix with n >= 1
%           A is irreducible s.t. all its off-diagonal entries are nonzero
%           tol is a finite positive number and ideally made small
% Time Complexity: O(n^2 * log_2(U - L / tol)), but better on AVERAGE
%                  Where n is number of rows/cols in the matrix A
%                  U and L are the upper and lower bounds on the
%                  eigenvalues from Gershgorin's circle theorem
%                  tol is the interval length tolerance
% Example Input: eigenTreeBisection([32,-16,0;-16,32,-16;0,-16,32], 1e-10)
% Example Output: [9.3726, 32.0000, 54.6274]

% Default arguments if none assigned
arguments
    A
    tol = 1e-5
end

% Use Gershgorin's Circle Theorem to find the bounds on all eigenvalues of A
circleCentres = diag(A);
b = abs(diag(A, -1));

% a_{1,1} circle has radius prop. to b_{1}, 
% While a_{2,2} has radius prop. b_{1} + b_{2} and so on.
circleRadius = [0; b] + [b; 0]; % This form guarantees the above works

% Take the minimums & maximums (on the real axis) of the union of all n-disks
lowerBound = min(circleCentres - circleRadius);
upperBound = max(circleCentres + circleRadius);

n = size(A, 1);
x = zeros(n, 1); % Initialize eigenvalue vector

% By Gershgorin's Theorem, all n eigenvalues of A lie within [lowerBound, upperBound]
x = splitInterval(A, lowerBound, upperBound, 0, n, tol, x);

end


function x = splitInterval(A, L, U, countL, countU, tol, x)
% x = splitInterval(A, L, U, countL, countU, tol, x) 
%   Returns the eigenvalues of A in increasing order, by recursively 
%   bisecting the search interval [L, U] until an interval length tolerance 
%   tol is reached.
%   countL and countU contain the number of eigenvalues of A below L and U
%   respectively.
% Requires: A is a symmetric tridiagonal square nxn matrix with n >= 1
%           A is irreducible s.t. all its off-diagonal entries are nonzero
%           L is a finite positive number
%           U is a finite positive number with U >= L
%           tol is a finite positive number and ideally made small
%           x is an array of size n

% Default arguments if none assigned
arguments
    A
    L
    U
    countL
    countU
    tol = 1e-5
    n
end

% How many eigenvalues in the current interval
numEigenvalues = countU - countL;

% If none, return
if numEigenvalues == 0
    return;
end
 

% Once interval tolerance is reached, stop dividing further
if (U - L) <= tol
    % Why assign countL+1:countU? It is because all remaining eigenvalues
    % therefore are below the interval tolerance and thus seen as
    % equivalent by the algorithm.
    x(countL+1:countU) = L + (U-L)/2; 
    return;
end

% Calculate midpoint after tolerance check to avoid catastrophic cancellation
midPoint = (U - L)/2;

% Otherwise, bisect the interval
newEndpoint = L + midPoint;

% If midPoint or L is below machine precision
if newEndpoint == L || newEndpoint == U
    error('Too small of a tolerance tol chosen');
end

% Count # of eigenvalues below the newly updated endpoint
countNew = sturmCounter(A, newEndpoint);

% Recursively work on the left interval
if countL < countNew % There are more eigenvalues in countNew
    x = splitInterval(A, L, newEndpoint, countL, countNew, tol, x); % set U = newEndpoint
end

% Recursively work on the right interval
if countU > countNew % There are less eigenvalues in countNew
    x = splitInterval(A, newEndpoint, U, countNew, countU, tol, x); % Set L = newEndpoint
end

end