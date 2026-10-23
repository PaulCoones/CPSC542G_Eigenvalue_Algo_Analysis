%{ 
    eigenJacobi.m provides a standalone function for computing the eigenvalues
    of a matrix A using the Jacobi eigenvalue algorithm (T.B. Chapter 30).
%}

function [D] = eigenJacobi(A)
    % x = eigenBisection(A) Returns the .....
    % Requires: A is a (real?) symmetric matrix
    % Time Complexity: ???
    % Output: A diagonalized matrix, where the diagonal elements correspond
    %         to the eigenvalues of A.
    % Example input: eigenJacobi([1, 2; 3, 7])
    % Example output: something?, compute it
    
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
        [j, k] = ind2sub(size(A), idx);
    
        % Formulate a Jacobi Matrix to zero out the largest off-diag element.
        
        % Calculate theta for jacobi matrix
        a_jj = A(j, j);
        a_kk = A(k, k);
        a_jk = A(j, k);

        if (a_jj - a_kk < tol)
            theta = pi / 4;
        else
            theta = 0.5*atan((2*a_jk)/(a_jj - a_kk));
        end

        % define Jacobi entries
        J_diag = cos(theta);
        J_jk = sin(theta);
        J_kj = -sin(theta);

        J = eye(n, n);
        J(j, j) = J_diag;
        J(k, k) = J_diag;
        J(j, k) = sin(theta);
        J(k, j) = -sin(theta);

        A = J'*A*J;

    end

    D = A
end