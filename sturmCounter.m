function signCount = sturmCounter(A, shift)
% signCount = sturmCounter(A, shift) Returns the number of sign changes
% in the Sturm sequence constructed from A with a provided shift.
% Requires: A is a symmetric tridiagonal nxn square matrix with n >= 1
%           A is irreducible s.t. all its off-diagonal entries are nonzero
%           shift is a number
% Time Complexity: O(n) s.t. n = size(A, 1)
% Space Complexity: O(n) s.t. n = size(A, 1)
% Example input: sturmCounter([32,-16,0;-16,32,-16;0,-16,32], 10)
% Example output: 1, which tells us the matrix A has 1 eigenvalue below 10

n = size(A, 1);
p = zeros(n+1, 1);     % We want sturm sequence from 0 to n = n+1 elements
p(1) = 1;              % Define p_0 = p(1) = 1
p(2) = A(1,1) - shift; % p_1 = a_{1,1} - alpha

for k = 2:n
    % Ref: Trefethen & Bau's Numerical Linear Algebra Equation (30.9)
    p(k+1) = ((A(k, k) - shift)*p(k)) - ((A(k, k-1)^2) * p(k-1));
end

% Count the number of sign changes in the shifted sturm sequence
signCount = sum((p(1:end-1) >= 0 & p(2:end) < 0) | ... 
    (p(1:end-1) <= 0 & p(2:end) > 0));

end