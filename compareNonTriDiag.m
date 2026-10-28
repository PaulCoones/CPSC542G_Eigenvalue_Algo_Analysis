%{ 
    compareNonTriDiag.m Is a program which compares exactly three robust
    (eigenvalue) algorithms against one another, by comparing their
    execution times and errors in their computed eigenvalues within a table.

    Hit run to run the program and display the table.
%}

% Clean workspace
clearvars;
close all;

% Test parameters (ensure length(nValues) = length(LValues))
nValues = [32];   % Add n tests here
LValues = [33]; % Add corresponding L values here

% Caclulate how many (n, L) pair tests we are performing
numTests = length(nValues);

% Initalize the three algorithms to compare
algorithmQRIter = @stdQRIteration; % QR Iteration
algorithmShiftQR = @qreig;         % Shifted QR Iteration
algorithmJacobi = @eigenJacobi;    % Jacobi

% Initialize errors
errorQRIter = zeros(numTests, 1);
errorShiftQR = zeros(numTests, 1);
errorJacobi = zeros(numTests, 1);

% Initialize timings
timeQRIter = zeros(numTests, 1);
timeShiftQR = zeros(numTests, 1);
timeJacobi = zeros(numTests, 1);

% Benchmark and record errors
for i = 1:numTests

    % Iterate over all testing pairs (n, L)
    n = nValues(i);
    L = LValues(i);

    % Construct symmetric matrix
    rng(77);
    A = rand(n);      % [0, 1] interval
    Z = (A + A') / 2; % Average entries for symmetry

    % Compute reference eigenvalues using MATLAB's built-in eig
    lambdaExact = sort(eig(Z));

    % Compute approximate eigenvalues using QR Iteration
    tic;
    lambdaQRIter = algorithmQRIter(Z);
    timeQRIter(i) = toc;

    % Compute approximate eigenvalues using Shifted QR Iteration
    tic;
    lambdaShiftQR = algorithmShiftQR(Z);
    timeShiftQR(i) = toc;

    % Compute approximate eigenvalues using Jacobi
    tic;
    lambdaJacobi = algorithmJacobi(Z);
    timeJacobi(i) = toc;

    % Incase eigenvalues not sorted, sort eigenvalues in ascending order
    % I.e. Jacobi does not sort them.
    lambdaQRIter = sort(lambdaQRIter(:));
    lambdaShiftQR = sort(lambdaShiftQR(:));
    lambdaJacobi = sort(lambdaJacobi(:));

    % Compute the errors in the computed eigenvalues and the reference
    % eigenvalues using the 2-norm
    errorQRIter(i) = norm(lambdaExact - lambdaQRIter, 2);
    errorShiftQR(i) = norm(lambdaExact - lambdaShiftQR, 2);
    errorJacobi(i) = norm(lambdaExact - lambdaJacobi, 2);

end

% Tabulate the results
results = table( ...
    nValues, LValues, ...
    errorQRIter, errorShiftQR, errorJacobi, ...
    timeQRIter, timeShiftQR, timeJacobi, ...
    'VariableNames', ...
    {'n', 'L', ...
     'Error2Norm_QRIter', 'Error2Norm_shiftQR', 'Error2Norm_Jacobi', ...
     'Time_QRIter(in s)', 'Time_shiftQR(in s)', 'Time_Jacobi(in s)'});

fprintf('\n Algorithm Table Comparison \n\n');

% Format n and L without scientific notation
results.n = compose("%.2f", nValues);
results.L = compose("%.2f", LValues);

% Otherwise use short exponential notation
format short e

% Display the table
disp(results);
