% x = triDiagRepeatedEigens(n, interval_size, L) Returns a symmetric matrix
% similar to the finite differenced laplacian returned by triDiag(), with
% zeros added at each interval point to fix temperature. Essentially
% subdivides the original heat rod problem into multiple sub-intervals so
% that eigenvalues repeat.
%
% Inputs:
%   - n : the number of points in each interval
%   - interval_count: number of intervals, defines the number of times eigenvalues repeat. 
%   - L: length of each interval
% Outputs:
%   - (n*interval_size x n*interval_size) symmetric matrix similar in 
%     structure to triDiag() but with 0's at each interval point. 
function [Z, lambdaExact] = triDiagRepeatedEigens(n, interval_count, L)
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
    lam1 = (4 / h^2) * sin(k * pi / (2 * (n + 1))).^2;
    % Repeat for each interval
    lambdaExact = sort(repmat(lam1, interval_size, 1));
end