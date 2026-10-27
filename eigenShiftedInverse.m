%{ 
    eigenShiftedInverse.m provides a standalone function for computing all
    the Eigenvectors of a matrix A given the eigenvalues, using the shifted
    inverse method.
%}


function eigenvectors = eigenShiftedInverse(A, eigenvalues, tol, max_iters, eps)
% x = eigenShiftedInverse(A, eigenvalues) Returns the eigenvectors of A in
% the order given by eigenvalues utilizing the shifted inverse method. For
% inverse iteration, we solve linear system (A - uI)x = v_{k-1} followed
% by v_k = x / norm(x) iteratively until v converges to the eigenvector of
% u.
% Requires: A is a full nxn matrix with n >= 1
%           tol is finite positive number that's ideally small
%           max_iters is an upper bound on the number of iterations run
%           before giving up.
% Time Complexity: O(2/3n^4) for general square matrices, O(n^2)
%                  for tridiagonal matrices
%           n is number of rows/columns in the matrix A
%           Assuming our eigenvalues are good, we only require ~1 iteration
%           to get the correct eigenvalues.
% Example input: A = [32,-16,0;-16,32,-16;0,-16,32];
%                d = eigenJacobi(A, 10, 1e-10);
%                eigenShiftedInverse(A, d)
% Example output: ans = 3×3
%                 -5.0000e-01	-7.0711e-01	-5.0000e-01
%                 -7.0711e-01	-3.5903e-15	7.0711e-01
%                 5.0000e-01	-7.0711e-01	5.0000e-01
arguments
    A;
    eigenvalues;
    tol = 1e-10;
    max_iters = 1000;
end
n = size(A, 1);
I = eye(n);

eigenvectors = zeros(numel(eigenvalues), n);
% if any eigenvalues didn't converge, then we assume we chose a bad vector
% and try again.
for i=1:2
    for i=1:numel(eigenvalues)
        % if eigenvector was solved in previous iter, skip.
        if (norm(eigenvectors(i, :)) > tol)
            continue;
        end

        u = eigenvalues(i);
        v_0 = rand(n, 1);
        v_i = v_0 / norm(v_0);
        B = A - u*I;
        for k=1:max_iters
            x = B \ v_i;
            v_i = x / norm(x);
            
            if (norm(A*v_i - u*v_i) < tol)
                v_i
                eigenvectors(i, :) = x / norm(x);
                break;
            end
        end
    end
end


end