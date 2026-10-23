%{ 
    bisection.m provides a standalone function for computing all eigenvalues
    of a matrix A using the Bisection eigenvalue algorithm with an interval
    decision tree.
    Reference: Trefethen & Bau's Numerical Linear Algebra Chapter 30
%}

function x = eigenTreeBisection(A, tol)
% x = eigenBisection(A) Returns the eigenvalues of A, provided (something
% about tol var) by counting the number of sign changes in its sturm
% sequence .....
% Requires: A is a symmetric tridiagonal square nxn matrix with n >= 1
%           A is irreducible s.t. all its off-diagonal entries are nonzero
%           tol is a positive number and ideally made small
% Time Complexity: O(n^???)
% Example input: eigenBisection([1, 2; 3, 7], 1e-10)
% Example output: something?, compute it

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

% By construction, all n eigenvalues of A lie within [lowerBound, upperBound]
x = splitInterval(A, lowerBound, upperBound, 0, n, tol, x);

end


function x = splitInterval(A, L, U, countL, countU, tol, x)
% The interval contains countU-countL eigenvalues

numEigenvalues = countU - countL;

% No eigenvalues in this interval
if numEigenvalues == 0
    return;
end

% Stops subdividing an interval once reached desired tolerance
if (U - L) <= tol
    x(countL+1:countU) = L + (U-L)/2;
    return;
end

% Otherwise, Bisect the interval
midPoint = L + (U-L)/2;

% Prevent floating-point stagnation
if midPoint == L || midPoint == U
    error('Requested tolerance cannot be resolved.');
end

% Count # of eigenvalues below the midPoint 
countMid = sturmCounter(A, midPoint);

% Recursively work on the left interval
if countL < countMid % There are more eigenvalues in countMid
    x = splitInterval(A, L, midPoint, countL, countMid, tol, x); % set U = MidPoint
end

% Recursively work on the right interval
if countU > countMid % There are less eigenvalues in countMid
    x = splitInterval(A, midPoint, U, countMid, countU, tol, x); % Set L = midPoint
end

end
