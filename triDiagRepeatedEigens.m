% x = triDiagRepeatedEigens(n, interval_size, L) Returns a symmetric matrix
% similar to the discretized laplacian returned by triDiag(), with
% zeros added at each chosen interval point to fix the temperature to zero.
% Ultimately, this subdivides the original problem into multiple 
% sub-intervals on which eigenvalues repeat.
% Inputs:
%   - n : the number of points in each interval
%   - interval_count: number of intervals, defines the number of times eigenvalues repeat. 
%   - L: length of each interval
% Outputs:
%   - Z: A symmetric matrix of size (n*interval_size x n*interval_size) 
%     similar in  structure to triDiag() but with 0's at each interval point
%   - repeatedLambdas: The exact eigenvalues of the matrix Z
function [Z, repeatedLambdas] = triDiagRepeatedEigens(n, interval_count, L)
    N = n * interval_count;
    h = L / (n + 1);

    % Construct full tridiagonal matrix
    Z = tridiag(N, h * (N+1));
    
    % Zero out off-diag elements at the interval joining points
    % Creates repeating tridiag structure
    for j = 1:interval_count-1
        last = j*n;
        Z(last, last+1) = 0;
        Z(last + 1, last) = 0;
    end

    % Exact eigenvalues of one segment
    k = (1:n)';
    originalLambdas = (4 / h^2) * sin(k * pi / (2 * (n + 1))).^2;

    % Repeat for each interval
    repeatedLambdas = sort(repmat(originalLambdas, interval_size, 1));
end