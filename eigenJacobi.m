%{ 
    eigenJacobi.m provides a standalone function for computing the eigenvalues
    of a matrix A using the Jacobi eigenvalue algorithm (T.B. Chapter 30).
%}

function [D] = eigenJacobi(A)
% x = eigenBisection(A) Returns the .....
% Requires: A is a (real?) symmetric matrix
% Time Complexity: O(kn^3) 
% Output: A diagonalized matrix, where the diagonal elements correspond
%         to the eigenvalues of A.
% Example input: eigenJacobi([1, 2; 2, 7])
% Example output: 

% Algorithm
[n, n] = size(A);

% Define a max number of iterations (k = 1000?)
max_iter = 1000;
tol = 1e-14; % floating point precision stops at 1e-16

% At each step:
for k=1:max_iter
    % Check if A is diagonal (stop if it is)
    off = norm(A - diag(diag(A)), "fro");

    if (off < tol)
        break;
    end

    % Determine the largest off-diag value of A

    % A is symmetric so we only need to consider the upper half
    B = abs(triu(A, 1)); % setting 1 discards the diagonal elements.
    [val, idx] = max(B(:));
    [p, q] = ind2sub(size(A), idx);

    % Formulate a Jacobi Matrix to zero out the largest off-diag element.

    % Calculate theta for jacobi matrix
    a_pp = A(p, p);
    a_qq = A(q, q);
    a_pq = A(p, q);

    if (a_pp - a_qq < tol)
        theta = pi / 4;
    else
        theta = 0.5*atan((2*a_pq)/(a_qq - a_pp));
    end

    % define Jacobi entries
    J_diag = cos(theta);
    J_pq = sin(theta);
    J_qp = -sin(theta);

    J = eye(n, n);
    J(p, p) = J_diag;
    J(q, q) = J_diag;
    J(p, q) = J_pq;
    J(q, p) = J_qp;

    A = J'*A*J;

end

D = A
end