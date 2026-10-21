%{ 
    bisection.m provides a standalone function for computing the eigenvalues
    of a matrix A using the Bisection eigenvalue algorithm (T.B. Chapter 30).
%}

function x = eigenBisection(A)
% x = eigenBisection(A) Returns the .....
% Requires: A is a (symmetric?) tridiagonal matrix
%           A is non-zero along any diagonal elements
% Time Complexity: O(n^2)
% Example input: eigenBisection([1, 2; 3, 7])
% Example output: something?, compute it