%{ 
    eigenParBisection.m provides a standalone function for computing all eigenvalues
    of a matrix A using the Bisection eigenvalue algorithm in naively and in Parallel.
    Reference: Trefethen & Bau's Numerical Linear Algebra Chapter 30
%}

function x = eigenParBisection(A, tol)
% x = eigenParBisection(A) Returns the eigenvalues of A in increasing order, 
%   by applying the Bisection eigenvalue algorithm for each eigenvalue 
%   until an interval length tolerance tol is reached.
%   This implementation is parallelized and naive, as for each eigenvalue,
%   we restart from the largest possible interval in which the eigenvalue
%   may be contained, however we can possibly do better, see eigenTreeBisection.m
% Requires: A is a real symmetric tridiagonal square nxn matrix with n >= 1
%           A is irreducible s.t. all its off-diagonal entries are nonzero
%           tol is a real and finite positive number, ideally made small
% Time Complexity: O(n^2/p * log_2(U - L / tol))
%                  Where n is number of rows/cols in the matrix A
%                  p is the number of workers
%                  U and L are the upper and lower bounds on the
%                  eigenvalues from Gershgorin's circle theorem
%                  tol is the interval length tolerance
% Example Input: eigenParBisection([32,-16,0;-16,32,-16;0,-16,32], 1e-10)
% Example Output: [9.3726, 32.0000, 54.6274]

% Default arguments if none assigned
arguments
    A
    tol = 1e-5
end

% Use Gershgorin's Circle Theorem to find the bounds on all eigenvalues of A
circleCentres = diag(A);
b = abs(diag(A, -1));

% a_{1,1} circle has radius prop. to b_{1}, 
% While a_{2,2} has radius prop. b_{1} + b_{2} and so on.
circleRadius = [0; b] + [b; 0]; % This form guarantees the above works

% Take the minimums & maximums (on the real axis) of the union of all n-disks
lowerBound = min(circleCentres - circleRadius);
upperBound = max(circleCentres + circleRadius);

n = size(A, 1);
x = zeros(n, 1); % Initialize eigenvalue vector

% Each eigenvalue is found independently, by bisecting the initial interval
% given by the bounds provided by Gershgorin's circle theorem.
% This results in a naive implementation, as we restart from the initial
% interval for each new eigenvalue we wish to compute.
% However, with this naivety, we can find each eigenvalue independently in parallel.
parfor i = 1:n
    L = lowerBound;
    U = upperBound;
    intervalLen = U - L;

    while intervalLen > tol

        midPoint = intervalLen / 2;
        curShift = L + midPoint;

        % If midPoint or L is below machine precision
        if curShift == L || curShift == U
            break; % Break out of the loop
        end

        eigenCount = sturmCounter(A, curShift); % O(n) complexity call

        % Bisection
        if eigenCount < i
            L = curShift;
        else
            U = curShift;
        end

        intervalLen = U - L; % Update new interval length

    end

    % Estimate this j-th eigenvalue upon hitting interval tol
    x(i) = L + (U - L)/2;

end

end

% %CONSOLE TEST
% if isempty(gcp('nocreate'))
%     parpool('Processes', 8);
% end
% 
% Z=triDiag(7000,1);
% 
% %Run parallel bisection
% tic;
% x1 = eigenParBisection(Z);
% elapsedTime1 = toc;
% fprintf('Parallel execution time: %.6f seconds\n', elapsedTime1);
% 
% tic;
% x2 = eigenBisection(Z);
% elapsedTime2 = toc;
% fprintf('Naive Sequential execution time: %.6f seconds\n', elapsedTime2);
% 
% tic;
% x3 = eigenTreeBisection(Z);
% elapsedTime3 = toc;
% fprintf('Good Sequential execution time: %.6f seconds\n', elapsedTime3);
% 
% tic;
% x4 = eigenJacobi(Z);
% elapsedTime4 = toc;
% fprintf('Jacobi execution time: %.6f seconds\n', elapsedTime4);
% delete(gcp('nocreate'))
