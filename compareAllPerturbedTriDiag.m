%{
    compareAllPerturbedTriDiag.m Is a program which compares all relevant
    (eigenvalue) algorithms against one another on a perturbed
    tridiagonal matrix, by comparing their execution times and 
    errors in their computed eigenvalues within a csv table.

    Hit run to run the program and save the results to a CSV table.
    Change the Test parameters and/or algorithm parameters to compare
    different (n, L) pairs and/or algorithms.

    Results are saved after every completed test.
    Local variables cleared at each iteration to save space.
%}

% Clean workspace
clearvars;
close all;

% Save benchmark results in a subfolder of the current directory
outputDir = fullfile(pwd, 'TabulatedBenchMarks');

% Create this directory if it doesnt yet exist
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Create a UNQIUE file name for the reuslts of this run according to the 
% current year+day+hour+minute+second.
fileName = ['compareAllPerturbedTriDiag' datestr(now, 'yyyymmdd_HHMMSS') '.csv'];

resultsFile = fullfile(outputDir, fileName);

% Test parameters (ensure length(nValues) = length(LValues))
nValues1 = 2.^(7);    % Test all six algorithms (n = 32 to 2048)
nValues2 = [2.^(2)];  % Test three algorithms (n = 4096 to 70000)

% Combine both sets of test parameters
nValues = [nValues1, nValues2]; % All n values to test
LValues = nValues + 1;          % Set h = L/(n+1) = 1 to normalize all tests

% Caclulate how many (n, L) pair tests we have to compute
numTests = length(nValues);

% Initalize the six algorithms to compare with anonymous functions
algorithmDivNConq = @eigenDivNConqNew;   % Divide and Conquer
algorithmTreeBi = @eigenTreeBisection;   % Tree Bisection
algorithmParBi = @eigenParBisection;     % Naive Parallel Bisection
algorithmNaiveBi = @eigenBisection;      % Naive Sequential Bisection
algorithmQRIter = @(A) stdQRIteration(A, 1e-8, 2500);   % QR Iteration
algorithmShiftQRIter = @(A) qreig(A, 1e-8, 2500, 0);    % Shifted QR Iteration
algorithmShiftQRIterRel = @(A) qreig(A, 1e-8, 2500, 1); % Shifted QR Iteration (Relative)
algorithmJacobi = @eigenJacobi;          % Jacobi

% Activate the parallel pool before performing any computations in parallel
if isempty(gcp('nocreate'))
    %parpool('Processes', 8); % Use cores
    parpool('Threads', 8);   % Use threads
end

fprintf('\nSaving the results to:\n %s \n', resultsFile);

% Benchmark and save results after every completed test
for i = 1:numTests

    % Iterate over all testing pairs (n, L)
    n = nValues(i);
    L = LValues(i);

    % Display progress of current test
    fprintf('\nRunning test %d of %d: n = %d, L = %.2f\n', ...
        i, numTests, n, L);

    % Perform the benchmark in a separate function so that its
    % large matrices and intermediate variables are released afterward
    resultsRow = runBenchmark(n, L, ...
        algorithmDivNConq, algorithmTreeBi, algorithmParBi, ...
        algorithmNaiveBi, algorithmQRIter, algorithmShiftQRIter, ...
        algorithmShiftQRIterRel, algorithmJacobi);

    % Create the CSV file for the first completed test
    if i == 1
        writetable(resultsRow, resultsFile);
    else
        % Add the completed test as a new row in the existing table
        writetable(resultsRow, resultsFile, 'WriteMode', 'append', ...
                   'WriteVariableNames', false);
    end

    % Confirmation that the current iteration results have been saved
    fprintf('Test %d of %d completed and saved successfully.\n', ...
        i, numTests);

    % Clear memory
    clear resultsRow;

end

% % Delete the parallel pool finishing all parallel executions
delete(gcp('nocreate'));

fprintf('\n All %d tests completed. \n', numTests);
fprintf('All results saved to: \n %s \n', resultsFile);



function CSVresultsRow = runBenchmark(n, L, ...
    algorithmDivNConq, algorithmTreeBi, algorithmParBi, ...
    algorithmNaiveBi, algorithmQRIter, algorithmShiftQRIter, ...
    algorithmShiftQRIterRel, algorithmJacobi)

    % Initialize all errors for the current test
    errorDivNConq = 0;
    errorTreeBi = 0;
    errorParBi = 0;
    errorNaiveBi = 0;
    errorQRIter = 0;
    errorShiftQRIter = 0;
    errorShiftQRIterRel = 0;
    errorJacobi = 0;

    % Initialize all timings for the current test
    timeDivNConq = 0;
    timeTreeBi = 0;
    timeParBi = 0;
    timeNaiveBi = 0;
    timeQRIter = 0;
    timeShiftQRIter = 0;
    timeShiftQRIterRel = 0;
    timeJacobi = 0;

    % Construct the tridiagonal matrix
    Z = triDiag(n, L);
    rng(7755);
    Z = Z + perturbTriDiag(Z);

    % Compute reference eigenvalues using MATLAB's eig solver
    lambdaExact = sort(eig(Z));

    % Compute the minimum spectral gap of the modified matrix
    minSpectralGap = min(diff(lambdaExact));

    % Compute approximate eigenvalues using Divide and Conquer
    [lambdaDivNConq, timeDivNConq] = timeAlgorithm(algorithmDivNConq, Z, n);

    % Incase eigenvalues not sorted, sort eigenvalues in ascending order
    lambdaDivNConq = sort(lambdaDivNConq(:));
    errorDivNConq = norm(lambdaExact - lambdaDivNConq, 2);
    clear lambdaDivNConq; % Clear memory


    % Compute approximate eigenvalues using Tree Bisection
    [lambdaTreeBi, timeTreeBi] = timeAlgorithm(algorithmTreeBi, Z, n);
    lambdaTreeBi = sort(lambdaTreeBi(:));
    errorTreeBi = norm(lambdaExact - lambdaTreeBi, 2);
    clear lambdaTreeBi;


    % Compute approximate eigenvalues using Naive Parallel Bisection
    [lambdaParBi, timeParBi] = timeAlgorithm(algorithmParBi, Z, n);
    lambdaParBi = sort(lambdaParBi(:));
    errorParBi = norm(lambdaExact - lambdaParBi, 2);
    clear lambdaParBi;


    % Only test all six algorithms for n = 32 to 2048
    if n <= 2048
        % Compute approximate eigenvalues using Naive Sequential Bisection
        [lambdaNaiveBi, timeNaiveBi] = timeAlgorithm(algorithmNaiveBi, Z, n);
        lambdaNaiveBi = sort(lambdaNaiveBi(:));
        errorNaiveBi = norm(lambdaExact - lambdaNaiveBi, 2);
        clear lambdaNaiveBi;


        % Compute approximate eigenvalues using QR Iteration
        [lambdaQRIter, timeQRIter] = timeAlgorithm(algorithmQRIter, Z, n);
        lambdaQRIter = sort(lambdaQRIter(:));
        errorQRIter = norm(lambdaExact - lambdaQRIter, 2);
        clear lambdaQRIter;


        % Compute approximate eigenvalues using Shifted QR Iteration
        [lambdaShiftQRIter, timeShiftQRIter] = timeAlgorithm( ...
                                            algorithmShiftQRIter, Z, n);
        lambdaShiftQRIter = sort(lambdaShiftQRIter(:));
        errorShiftQRIter = norm(lambdaExact - lambdaShiftQRIter, 2);
        clear lambdaShiftQRIter;


        % Compute approximate eigenvalues using Shifted QR Iteration (Relative)
        [lambdaShiftQRIterRel, timeShiftQRIterRel] = timeAlgorithm( ...
                                            algorithmShiftQRIterRel, Z, n);
        lambdaShiftQRIterRel = sort(lambdaShiftQRIterRel(:));
        errorShiftQRIterRel = norm(lambdaExact - lambdaShiftQRIterRel, 2);
        clear lambdaShiftQRIterRel;


        % Compute approximate eigenvalues using Jacobi
        [lambdaJacobi, timeJacobi] = timeAlgorithm(algorithmJacobi, Z, n);
        lambdaJacobi = sort(lambdaJacobi(:));
        errorJacobi = norm(lambdaExact - lambdaJacobi, 2);
        clear lambdaJacobi;
    end

    % Tabulate the results for the currents test
    CSVresultsRow = table( ...
        n, L, minSpectralGap, errorNaiveBi, errorParBi, errorTreeBi, ...
        errorQRIter, errorShiftQRIter, errorShiftQRIterRel, ...
        errorJacobi, errorDivNConq, ...
        timeNaiveBi, timeParBi, timeTreeBi, ...
        timeQRIter, timeShiftQRIter, timeShiftQRIterRel, ...
        timeJacobi, timeDivNConq, ...
        'VariableNames', ...
        {'n', 'L', 'MinSpectralGap', ...
         'Error2Norm_NaiveBisect', 'Error2Norm_ParBisect', ...
         'Error2Norm_TreeBisect', 'Error2Norm_QRIter', ...
         'Error2Norm_ShiftQRIter', 'Error2Norm_ShiftQRIterRel', ...
         'Error2Norm_Jacobi', 'Error2Norm_DivNConq', ...
         'Time_NaiveBisect(in s)', 'Time_ParBisect(in s)', ...
         'Time_TreeBisect(in s)', 'Time_QRIter(in s)', ...
         'Time_ShiftQRIter(in s)', 'Time_ShiftQRIterRel(in s)', ...
         'Time_Jacobi(in s)', 'Time_DivNConq(in s)'});

    % Clear intermediate vars for more space
    clear Z lambdaExact minSpectralGap k h;

end


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