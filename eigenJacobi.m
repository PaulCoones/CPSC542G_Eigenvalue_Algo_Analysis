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

max_sweeps = 100;
tol = 1e-15; % floating point precision stops at 1e-16

% Mentioned in TB Chapter 30. To avoid O(n^2) search, we iterate
% over sweeps of the the off-diag elements.
for kt=1:max_sweeps
    for j=1:n
        for k=j:n
            % Check if A is diagonal (stop if it is)
            off = norm(A - diag(diag(A)), "fro");

            if (off < tol)
                D = A;
                return;
            end

            % Calculate theta for jacobi matrix
            a_jj = A(j, j);
            a_kk = A(k, k);
            a_jk = A(j, k);
        
            if (a_jj - a_kk < tol)
                theta = pi / 4;
            else
                theta = 0.5*atan((2*a_jk)/(a_kk - a_jj));
            end
        
            % define Jacobi entries
            J_diag = cos(theta);
            J_jk = sin(theta);
            J_kj = -sin(theta);
        
            J = eye(n, n);
            J(j, j) = J_diag;
            J(k, k) = J_diag;
            J(j, k) = J_jk;
            J(k, j) = J_kj;
        
            A = J'*A*J;
        end
    end

end
warning("Jacobi eigenvalue method was not able to converge for matrx")
D = A;
end