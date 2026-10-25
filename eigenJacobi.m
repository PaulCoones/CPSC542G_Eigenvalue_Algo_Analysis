%{ 
    eigenJacobi.m provides a standalone function for computing all eigenvalues
    of a matrix A using the Jacobi eigenvalue algorithm (T.B. Chapter 30).
%}



function D = eigenJacobi(A, max_sweeps, tol)
% x = eigenJacobi(A) Returns the eigenvalues of A in increasing order
% utilizing the Jacobi eigenvalue method. Starting from the original matrix
% A we sequentially transform the rows and columns of A using Jacobi
% rotation matrices that 0 out the off-diagonal elements. After a
% sufficient number of applications the matrix converges to a diagonal
% matrix containing the eigenvalues of A
% Requires: A is a full symmetric square nxn matrix with n >= 1
%           max_sweeps is a non-zero positive integer that's ideally
%           greater than 1
%           tol is finite positive number that's ideally small
% Time Complexity: O(max_sweeps * n^3) 
%           max_sweeps is the upper bound on the number of iterations
%           in the algorithm.
%           n is number of rows/columns in the matrix A
% Example input: eigenJacobi([32,-16,0;-16,32,-16;0,-16,32], 10, 1e-10)
% Example output: [9.3726; 32.0000; 54.6274]

arguments
    A 
    max_sweeps=100
    tol=1e-5
end

n = size(A, 1);

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

            if (a_jk == 0)
                continue; % entry already 0
            end
        
            % % let t = tan(theta), x = tan(2*theta)
            % % x = (2t)/(1-t^2) <- double angle formula
            % % x - xt^2 - 2t = 0
            % % t^2 + (2/x)t - 1 = 0
            % 
            % % 1/x = 1/tan(2*\theta) = cot(2*\theta) = r
            % 
            % % t^2 + 2rt - 1 = 0
            % % (t + r)^2 - r^2 - 1 = 0
            % % (t + r)^2 = 1 + r^2
            % % t = -r +/- sqrt(1 + r^2)
            % 
            % % two solutions, so pick smallest
            % % if r >= 0 -> t = -r + sqrt(1 + r^2)
            % % if r < 0 -> t = -r - sqrt(1 + r^2)
            % 
            % % this version of t is numerically unstable
            % % t = (-r +/- sqrt(1 + r^2)) * (r +/- sqrt(1 + r^2))/(r +/- sqrt(1 + r^2))
            % % t = (-r^2 + 1 + r^2)/(r +/- sqrt(1 + r^2))
            % % t = 1/(r +/- sqrt(1+r^2))
            % 
            % % 1 + t^2 = 1/c^2 (c = cos(theta))
            % % c = 1/sqrt(1-t^2)
            % % s = tc (s = sin(theta))
            % 
            % r = (a_kk - a_jj) / (2*a_jk);
            % if r >= 0
            %     t = 1 / (r + sqrt(1 + r^2));
            % else
            %     t = 1 / (r - sqrt(1 + r^2));
            % end
            % 
            % % t != 0 (try and solve above formulas for t=0 if unsure)
            % c = 1 / sqrt(1 + t^2);
            % s = t * c;
            
            theta = 0.5*acot((a_kk - a_jj) / (2*a_jk));
            c = cos(theta);
            s = sin(theta);

            % define Jacobi entries
            J_diag = c;
            J_jk = s;
            J_kj = -J_jk;
        
            J_sub = [J_diag, J_jk; J_kj, J_diag];

            % only perform matrix computation on the target rows/columns
            A(:, [j, k]) = A(:, [j, k])*J_sub;
            A([j, k], :) = J_sub'*A([j, k], :);
            A(j, k) = 0;
            A(k, j) = 0;

        end
    end
end
D = sort(diag(A)); % To match the output of the Bisection method and others
end