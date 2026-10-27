%{ 
    stdQRIteration.m provides a standalone function for computing all 
    eigenvalues of a matrix A using QR Iteration (without shifts).
    Adapted from Dr. Chen Greifs qreig.m file for CPSC542G 2026WT1.
%}


function lambda = stdQRIteration(A, tol, max_iter)
% x = stdQRIteration(A) Returns the eigenvalues of A in increasing order, 
%   by applying QR Iteration until an tolerance tol is reached or a 
%   max number of iterations max_iter is attained.
% Requires: A is a square nxn matrix
%           max_iter is an integer, ideally made large >= 1e+4
%           tol is a finite positive number and ideally made small <= 1e-4
% Time Complexity: O(Mn^3) s.t. M = min(max_iter, i_{conv}) 
%                          where i is the number of iterations until the
%                          relative convergence criterion is met.
% Example Input: stdQRIteration(triDiag(3, 1), 1e-5, 10000)
% Example Output: [9.3726; 32.0000; 54.6274]

arguments
    A
    tol = 1e-5
    max_iter = 10000
end

% Compute a (relative) convergence criterion
convCriterion = tol * norm(A, 'fro'); 

i = 0;
% Perform QR iterations on the entire matrix
while i < max_iter

    % Check if A is close enough to being diagonal
    % By checking the worst case, hence Frobenius norm
    offDiag = A - diag(diag(A));
    if norm(offDiag, 'fro') <= convCriterion
        break; % If it is, break out and return
    end

    % Otherwise, compute the QR factorization of A
    [Q,R] = qr(A);

    % Apply the QR in reverse
    A = R*Q;

    % Increment
    i = i + 1;

end

% Return and sort the eigenvalues in increasing order
lambda = sort(diag(A));

end
