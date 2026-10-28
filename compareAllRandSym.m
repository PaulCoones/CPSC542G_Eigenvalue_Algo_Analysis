%{
    compareAllRandSym.m Is a program which compares all relevant
    (eigenvalue) algorithms on a random symmetric matrix constructed from
    By comparing their execution times and errors in their computed 
    eigenvalues within a csv table.

    Hit run to run the program and save the results to a CSV table.
    Change the Test parameters and/or algorithm parameters to compare
    different matrix dimensions and/or algorithms.

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
fileName = ['directNaiveCompareAll_' datestr(now, 'yyyymmdd_HHMMSS') '.csv'];

resultsFile = fullfile(outputDir, fileName);

% Test parameters
nValues1 = 2.^(5:11);           % Test all seven algorithms (n = 32 to 2048)
nValues2 = [2.^(12:16), 70000]; % Test three algorithms (n = 4096 to 70000)

% Combine both sets of test parameters
nValues = [nValues1, nValues2]; % All n values to test

% Caclulate how many tests we have to compute
numTests = length(nValues);

% Initalize the seven algorithms to compare with anonymous functions
algorithmDivNConq = @eigenDivNConqNew;   % Divide and Conquer
algorithmTreeBi = @eigenTreeBisection;   % Tree Bisection
algorithmParBi = @eigenParBisection;     % Naive Parallel Bisection
algorithmNaiveBi = @eigenBisection;      % Naive Sequential Bisection
algorithmQRIter = @stdQRIteration;       % QR Iteration
algorithmShiftQRIter = @qreig;           % Shifted QR Iteration
algorithmJacobi = @eigenJacobi;          % Jacobi

% Activate the parallel pool before performing any computations in parallel
if isempty(gcp('nocreate'))
    %parpool('Processes', 8); % Use cores
    parpool('Threads', 8);   % Use threads
end

fprintf('\nSaving the results to:\n %s \n', resultsFile);

% Benchmark and save results after every completed test
for i = 1:numTests

    % Iterate over all testing dimensions
    n = nValues(i);

    % Display progress of current test
    fprintf('\nRunning test %d of %d: n = %d\n', ...
        i, numTests, n);

    % Perform the benchmark in a separate function so that its
    % large matrices and intermediate variables are released afterward
    resultsRow = runBenchmark(n, ...
        algorithmDivNConq, algorithmTreeBi, algorithmParBi, ...
        algorithmNaiveBi, algorithmQRIter, algorithmShiftQRIter, ...
        algorithmJacobi);

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


function CSVresultsRow = runBenchmark(n, algorithmDivNConq, algorithmTreeBi, ...
    algorithmParBi, algorithmNaiveBi, algorithmQRIter, algorithmShiftQRIter, ...
    algorithmJacobi)

    % Initialize all errors for the current test
    errorDivNConq = 0;
    errorTreeBi = 0;
    errorParBi = 0;
    errorNaiveBi = 0;
    errorQRIter = 0;
    errorShiftQRIter = 0;
    errorJacobi = 0;

    % Initialize all timings for the current test
    timeDivNConq = 0;
    timeTreeBi = 0;
    timeParBi = 0;
    timeNaiveBi = 0;
    timeQRIter = 0;
    timeShiftQRIter = 0;
    timeJacobi = 0;

    % Construct the random symmetric matrix
    rng(77);
    A = rand(n);
    Z = (A + A') / 2; % Average entries for symmetry.
    clear A;

    % Compute all reference eigenvalues using MATLAB eig
    lambdaExact = sort(eig(Z));

    % Compute the symmetric tridiagonal matrix using Householder
    [T, timeHouseEig] = timeAlgorithm(@houseeig, Z, n);

    % Test Divide and Conquer, Tree Bisection and Parallel Bisection
    % for all values of n!
    % Compute approximate eigenvalues using Divide and Conquer
    [lambdaDivNConq, timeDivNConq] = timeAlgorithm(algorithmDivNConq, T, n);
    timeDivNConq = timeDivNConq + timeHouseEig;
    % Incase eigenvalues not sorted, sort eigenvalues in ascending order
    lambdaDivNConq = sort(lambdaDivNConq(:));
    errorDivNConq = norm(lambdaExact - lambdaDivNConq, 2);
    clear lambdaDivNConq; % Clear memory

    % Compute approximate eigenvalues using Tree Bisection
    [lambdaTreeBi, timeTreeBi] = timeAlgorithm(algorithmTreeBi, T, n);
    timeTreeBi = timeTreeBi + timeHouseEig;
    lambdaTreeBi = sort(lambdaTreeBi(:));
    errorTreeBi = norm(lambdaExact - lambdaTreeBi, 2);
    clear lambdaTreeBi;

    % Compute approximate eigenvalues using Naive Parallel Bisection
    [lambdaParBi, timeParBi] = timeAlgorithm(algorithmParBi, T, n);
    timeParBi = timeParBi + timeHouseEig;
    lambdaParBi = sort(lambdaParBi(:));
    errorParBi = norm(lambdaExact - lambdaParBi, 2);
    clear lambdaParBi;

    % Test all seven algorithms for n = 32 to 2048
    if n <= 2048
        % Compute approximate eigenvalues using Naive Sequential Bisection
        [lambdaNaiveBi, timeNaiveBi] = timeAlgorithm(algorithmNaiveBi, T, n);
        timeNaiveBi = timeNaiveBi + timeHouseEig;
        lambdaNaiveBi = sort(lambdaNaiveBi(:));
        errorNaiveBi = norm(lambdaExact - lambdaNaiveBi, 2);
        clear lambdaNaiveBi;

        % Compute approximate eigenvalues using QR Iteration
        [lambdaQRIter, timeQRIter] = timeAlgorithm(algorithmQRIter, T, n);
        timeQRIter = timeQRIter + timeHouseEig;
        lambdaQRIter = sort(lambdaQRIter(:));
        errorQRIter = norm(lambdaExact - lambdaQRIter, 2);
        clear lambdaQRIter;

        % Compute approximate eigenvalues using Shifted QR Iteration
        [lambdaShiftQRIter, timeShiftQRIter] = timeAlgorithm( ...
                                            algorithmShiftQRIter, T, n);
        timeShiftQRIter = timeShiftQRIter + timeHouseEig;
        lambdaShiftQRIter = sort(lambdaShiftQRIter(:));
        errorShiftQRIter = norm(lambdaExact - lambdaShiftQRIter, 2);
        clear lambdaShiftQRIter;

        % Compute approximate eigenvalues using Jacobi
        [lambdaJacobi, timeJacobi] = timeAlgorithm(algorithmJacobi, Z, n);
        lambdaJacobi = sort(lambdaJacobi(:));
        errorJacobi = norm(lambdaExact - lambdaJacobi, 2);
        clear lambdaJacobi;
    end

    % Tabulate the results for the currents test
    CSVresultsRow = table( ...
        n, errorNaiveBi, errorParBi, errorTreeBi, ...
        errorQRIter, errorShiftQRIter, errorJacobi, errorDivNConq, ...
        timeNaiveBi, timeParBi, timeTreeBi, ...
        timeQRIter, timeShiftQRIter, timeJacobi, timeDivNConq, ...
        'VariableNames', ...
        {'n', 'Error2Norm_NaiveBisect', 'Error2Norm_ParBisect', ...
         'Error2Norm_TreeBisect', 'Error2Norm_QRIter', ...
         'Error2Norm_ShiftQRIter', 'Error2Norm_Jacobi', ...
         'Error2Norm_DivNConq', ...
         'Time_NaiveBisect(in s)', 'Time_ParBisect(in s)', ...
         'Time_TreeBisect(in s)', 'Time_QRIter(in s)', ...
         'Time_ShiftQRIter(in s)', 'Time_Jacobi(in s)', 'Time_DivNConq(in s)'});

    % Clear intermediate vars for more space
    clear A Z T lambdaExact timeHouseEig;

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
