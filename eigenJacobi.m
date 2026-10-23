%{ 
    eigenJacobi.m provides a standalone function for computing the eigenvalues
    of a matrix A using the Jacobi eigenvalue algorithm (T.B. Chapter 30).
%}

function D = eigenJacobi(A, max_sweeps, tol)
% x = eigenBisection(A) Returns the .....
% Requires: A is a real symmetric matrix
% Time Complexity: O(max_sweeps * n^3) 
% Output: A diagonalized matrix, where the diagonal elements correspond
%         to the eigenvalues of A.
% Example input: ?? ADD ME ??
% Example output: ?? ADD ME ??

arguments
    A 
    max_sweeps=100
    tol=1e-5
end

% We can assume A is not sparse? Requires statement
if issparse(A)
    A = full(A);
end

% Algorithm
n = size(A, 1);
%J = eye(n, n);

mask = true(1, size(A, 1));

% Mentioned in TB Chapter 30. To avoid an O(n^2) search, we iterate
% over sweeps of the the off-diag elements.
for kt=1:max_sweeps
    % Check if A is diagonal and stop if it is
    off = norm(A - diag(diag(A)), "fro");
    if (off < tol)
        break;
    end
    for j=1:n-1
        for k=j+1:n

            % Calculate theta for jacobi matrix
            a_jj = A(j, j);
            a_kk = A(k, k);
            a_jk = A(j, k);
        
            if (a_jj == a_kk)
                theta = pi / 4;
            else
                theta = 0.5*atan((2*a_jk)/(a_kk - a_jj));
            end
        
            % define Jacobi entries
            J_diag = cos(theta);
            J_jk = sin(theta);
            J_kj = -J_jk;
        
            %J(j, j) = J_diag;
            %J(k, k) = J_diag;
            %J(j, k) = J_jk;
            %J(k, j) = J_kj;
            J_sub = [J_diag, J_jk; J_kj, J_diag];

            % only perform matrix computation on the target rows/columns
            A(:, [j, k]) = A(:, [j, k])*J_sub;
            A([j, k], :) = J_sub'*A([j, k], :);
            A(j, k) = 0;
            A(k, j) = 0;
            %A = J'*A*J;

            % zero out elements so J goes back to identity
            %J(j, j) = 1;
            %J(k, k) = 1;
            %J(j, k) = 0;
            %   J(k, j) = 0;
        end
    end
end
D = diag(A); % To match the output of the Bisection method and others
end