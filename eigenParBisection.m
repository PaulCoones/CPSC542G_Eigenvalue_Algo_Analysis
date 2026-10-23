%{ 
    bisection.m provides a standalone function for computing all eigenvalues
    of a matrix A using the Bisection eigenvalue algorithm in Parallel.
    Reference: Trefethen & Bau's Numerical Linear Algebra Chapter 30
%}

function x = eigenParBisection(A, tol)
% x = eigenBisection(A) Returns the eigenvalues of A, provided (something
% about tol var) by counting the number of sign changes in its sturm
% sequence .....
% Requires: A is a symmetric tridiagonal square nxn matrix with n >= 1
%           A is irreducible s.t. all its off-diagonal entries are nonzero
%           tol is a positive number and ideally made small
% Time Complexity: O(n^???)
% Example input: eigenBisection([1, 2; 3, 7], 1e-10)
% Example output: something?, compute it

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
parfor i = 1:n
    L = lowerBound;
    U = upperBound;
    intervalLen = U - L;

    while intervalLen > tol

        midPoint = intervalLen / 2;
        curShift = L + midPoint;

        % Incase midPoint is below machine Precision and we do not move,
        % or L is negligible compared to U.
        if curShift == L || curShift == U
            break;
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

% CONSOLE TEST
% if isempty(gcp('nocreate'))
%     parpool('Processes', 32);
% end
% 
% Z=triDiag(7000,1);
% 
% % Run parallel bisection
% tic;
% x1 = eigenParBisection(Z, 1e-5);
% elapsedTime1 = toc;
% fprintf('Parallel execution time: %.6f seconds\n', elapsedTime1);
% 
% tic;
% x2 = eigenBisection(Z, 1e-5);
% elapsedTime2 = toc;
% fprintf('Naive Sequential execution time: %.6f seconds\n', elapsedTime2);
% 
% tic;
% x3 = eigenTreeBisection(Z, 1e-5);
% elapsedTime3 = toc;
% fprintf('Good Sequential execution time: %.6f seconds\n', elapsedTime3);