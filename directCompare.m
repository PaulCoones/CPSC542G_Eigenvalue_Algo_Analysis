% FIX ME PAUL
% FIX ME PAUL

% Use to directly compare two algorithms on specific sets of
% n and L values.

% Clean workspace
clear all;
close all;

% =====================================================
% INITIALIZATION
% =====================================================

% Test parameters (ensure length(nValues) = length(LValues)
nValues = [25; 50; 75; 100; 500]; % Add n tests here
LValues = [1; 1; 1; 1; 501];      % Add corresponding L values here
% For Bisection, if nValue = n, LValue = n+1 enables convergence

numTests = length(nValues);

% Algorithms to compare 
algorithmX = @eigenBisection; % Change X to, for example, eigenParBisection or other
algorithmY = @eigenJacobi;    % Change Y to, for example, eigenTreeBisection or other

% Initialize error vectors
errorX = zeros(numTests, 1);
errorY = zeros(numTests, 1);

% Initialize execution time vectors
timeX = zeros(numTests, 1);
timeY = zeros(numTests, 1);

% =====================================================
% COMPUTE EIGENVALUES AND ERRORS
% =====================================================

for i = 1:numTests

    n = nValues(i);
    L = LValues(i);

    % Construct tridiagonal matrix
    Z = triDiag(n, L);

    % Compute exact analytical eigenvalues
    h = L/(n+1);
    k = (1:n)';

    lambdaExact = (4/h^2) * ...
        sin(k*pi/(2*(n+1))).^2;

    % Compute approximate eigenvalues using Algorithm X
    tic;
    lambdaX = algorithmX(Z);
    timeX(i) = toc;

    % Compute approximate eigenvalues using Algorithm Y
    tic;
    lambdaY = algorithmY(Z);
    timeY(i) = toc;

    % Incase eigenvalues not sorted, sort eigenvalues in ascending order
    lambdaX = sort(lambdaX(:));
    lambdaY = sort(lambdaY(:));

    % Compute 2-norm errors
    errorX(i) = norm(lambdaExact - lambdaX, 2);
    errorY(i) = norm(lambdaExact - lambdaY, 2);

end

% =====================================================
% DISPLAY RESULTS
% =====================================================

results = table( ...
    nValues, LValues, ...
    errorX, errorY, timeX, timeY, ...
    'VariableNames', ...
    {'n', 'L', ...
     'Error2Norm_X', 'Error2Norm_Y', ...
     'Time_X_s', 'Time_Y_s'});

fprintf('\n=== EIGENVALUE ALGORITHM COMPARISON ===\n\n');

% Format n and L without scientific notation
results.n = compose("%.2f", nValues);
results.L = compose("%.2f", LValues);

format short e
disp(results);