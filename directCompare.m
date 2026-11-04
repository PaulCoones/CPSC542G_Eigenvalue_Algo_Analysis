%{ 
    directCompare.m Is a program which compares exactly two (eigenvalue)
    algorithms against one another, by comparing their execution times
    and errors in their computed eigenvalues within a table.

    Hit run to run the program and display the table.
    Change the Test parameters and/or algorithmX and/or algorithmY
    parameters to compare different (n, L) pairs and/or algorithms.
%}

% Clean workspace
clearvars;
close all;

% Test parameters (ensure length(nValues) = length(LValues)
nValues = [32; 64; 64; 64; 64; 128; 256]; % Add n tests here
LValues = [33; 0.01; 1; 65; 129; 129; 257];      % Add corresponding L values here
% For Bisection, if nValue = n, LValue = n+1 enables error convergence
% with certainty. As otherwise, values in the triDiag matrix grow
% without bound with increasing n.

% Caclulate how many (n, L) pair tests we are performing
numTests = length(nValues);

% Initalize the two algorithms to compare
algorithmX = @(A) eigenJacobi(A, 100, 1e-5, 1); % Change X to, for example, eigenParBisection or other
algorithmY = @(A) eigenJacobi(A, 100, 1e-5, 0);    % Change Y to, for example, eigenTreeBisection or other

% Initialize errors
errorX = zeros(numTests, 1);
errorY = zeros(numTests, 1);

% Initialize timings
timeX = zeros(numTests, 1);
timeY = zeros(numTests, 1);

% % Activate the parallel pool before performing computations in parallel
%parpool('Processes', 8); % Use cores
parpool('Threads', 8);    % Use threads

% Benchmark and record errors
for i = 1:numTests

    % Iterate over all testing pairs (n, L)
    n = nValues(i);
    L = LValues(i);

    % Construct tridiagonal matrix
    Z = triDiag(n, L);

    % Compute exact eigenvalues from the known analytical solution
    h = L/(n+1);
    k = (1:n)';
    lambdaExact = (4/h^2) * ...
        sin(k*pi/(2*(n+1))).^2;

    % Compute approximate eigenvalues using Algorithm X
    %tic;
    lambdaX = algorithmX(Z);
    [lambdaX, xTime] = timeAlgorithm(algorithmX, Z, )
    %timeX(i) = toc;

    % Compute approximate eigenvalues using Algorithm Y
    %tic;
    lambdaY = algorithmY(Z);
    %timeY(i) = toc;

    % Incase eigenvalues not sorted, sort eigenvalues in ascending order
    % I.e. Jacobi does not sort them.
    lambdaX = sort(lambdaX(:));
    lambdaY = sort(lambdaY(:));

    % Compute the errors in the computed eigenvalues and the exact
    % eigenvalues using the 2-norm
    errorX(i) = norm(lambdaExact - lambdaX, 2);
    errorY(i) = norm(lambdaExact - lambdaY, 2);

end

% % Delete the parallel pool after having finished parallel executions
delete(gcp('nocreate'));

% Tabulate the results
results = table( ...
    nValues, LValues, ...
    errorX, errorY, timeX, timeY, ...
    'VariableNames', ...
    {'n', 'L', ...
     'Error2Norm_AlgoX', 'Error2Norm_AlgoY', ...
     'Time_AlgoX(in s)', 'Time_AlgoY(in s)'});

fprintf('\n Algorithm Table Comparison \n\n');

% Format n and L without scientific notation
results.n = compose("%.2f", nValues);
results.L = compose("%.2f", LValues);

% Otherwise use short exponential notation
format short e

% Display the table
disp(results);

function [eigens, measuredTime] = timeAlgorithm(algorithm, Z, n)
if n < 4096 % Use timeit for n < 4096 since runtimes are relatively small
    measuredTime = timeit(@() algorithm(Z));
    eigens = algorithm(Z);
else % Use tic/toc for n >= 4096 since runtimes become large
    tic;
    eigens = algorithm(Z);
    measuredTime = toc;
end
end