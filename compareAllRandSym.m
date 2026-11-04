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
fileName = ['finalLARGENCompareAllRandSym_' datestr(now, 'yyyymmdd_HHMMSS') '.csv'];

resultsFile = fullfile(outputDir, fileName);

% Test parameters
nValues1 = 2.^(2);           % Test all seven algorithms (n = 32 to 2048)
nValues2 = 2.^(11:16); % Test three algorithms (n = 4096 to 65536)

% Combine both sets of test parameters
nValues = [nValues1, nValues2]; % All n values to test

% Caclulate how many tests we have to compute
numTests = length(nValues);

% Initalize the seven algorithms to compare with anonymous functions
algorithmDivNConq = @eigenDivNConqNew;   % Divide and Conquer
algorithmTreeBi = @eigenTreeBisection;   % Tree Bisection
algorithmParBi = @eigenParBisection;     % Naive Parallel Bisection
algorithmNaiveBi = @eigenBisection;      % Naive Sequential Bisection
algorithmQRIter = @(A) stdQRIteration(A, 1e-8, 2500);       % QR Iteration
algorithmShiftQRIter = @(A) qreig(A, 1e-8, 2500, 0); % Shifted QR Iteration
algorithmShiftQRIterRel = @(A) qreig(A, 1e-8, 2500, 1); % Shifted QR Iteration w/ rel tol 
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


function CSVresultsRow = runBenchmark(n, algorithmDivNConq, algorithmTreeBi, ...
    algorithmParBi, algorithmNaiveBi, algorithmQRIter, algorithmShiftQRIter, ...
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

    % Construct the random symmetric matrix
    rng(77);
    A = rand(n);
    Z = (A + A') / 2; % Average entries for symmetry.
    clear A;

    % Compute all reference eigenvalues using MATLAB eig
    lambdaExact = sort(eig(Z));

    % Compute the minimum spectral gap from the reference eigenvalues
    minSpectralGap = min(diff(lambdaExact));

    % 2-norm of the reference eigenvalues (used for the relative errors)
    normExact = norm(lambdaExact, 2);

    % Compute the symmetric tridiagonal matrix using MATLABs hess
    [T, timeTridiagonalize] = timeAlgorithm(@hess, Z, n);

    % Test Divide and Conquer, Tree Bisection and Parallel Bisection
    % for all values of n!
    % Compute approximate eigenvalues using Divide and Conquer
    [lambdaDivNConq, timeDivNConq] = timeAlgorithm(algorithmDivNConq, T, n);
    timeDivNConq = timeDivNConq + timeTridiagonalize;
    % Incase eigenvalues not sorted, sort eigenvalues in ascending order
    lambdaDivNConq = sort(lambdaDivNConq(:));
    errorDivNConq = norm(lambdaExact - lambdaDivNConq, 2);
    clear lambdaDivNConq; % Clear memory

    % Compute approximate eigenvalues using Tree Bisection
    [lambdaTreeBi, timeTreeBi] = timeAlgorithm(algorithmTreeBi, T, n);
    timeTreeBi = timeTreeBi + timeTridiagonalize;
    lambdaTreeBi = sort(lambdaTreeBi(:));
    errorTreeBi = norm(lambdaExact - lambdaTreeBi, 2);
    clear lambdaTreeBi;

    % Compute approximate eigenvalues using Naive Parallel Bisection
    [lambdaParBi, timeParBi] = timeAlgorithm(algorithmParBi, T, n);
    timeParBi = timeParBi + timeTridiagonalize;
    lambdaParBi = sort(lambdaParBi(:));
    errorParBi = norm(lambdaExact - lambdaParBi, 2);
    clear lambdaParBi;

    % Test all seven algorithms for n = 32 to 2048
    if n <= 256
        % Compute approximate eigenvalues using Naive Sequential Bisection
        [lambdaNaiveBi, timeNaiveBi] = timeAlgorithm(algorithmNaiveBi, T, n);
        timeNaiveBi = timeNaiveBi + timeTridiagonalize;
        lambdaNaiveBi = sort(lambdaNaiveBi(:));
        errorNaiveBi = norm(lambdaExact - lambdaNaiveBi, 2);
        clear lambdaNaiveBi;

        % Compute approximate eigenvalues using QR Iteration
        [lambdaQRIter, timeQRIter] = timeAlgorithm(algorithmQRIter, T, n);
        timeQRIter = timeQRIter + timeTridiagonalize;
        lambdaQRIter = sort(lambdaQRIter(:));
        errorQRIter = norm(lambdaExact - lambdaQRIter, 2);
        clear lambdaQRIter;

        % Compute approximate eigenvalues using Shifted QR Iteration
        [lambdaShiftQRIter, timeShiftQRIter] = timeAlgorithm( ...
                                            algorithmShiftQRIter, T, n);
        timeShiftQRIter = timeShiftQRIter + timeTridiagonalize;
        lambdaShiftQRIter = sort(lambdaShiftQRIter(:));
        errorShiftQRIter = norm(lambdaExact - lambdaShiftQRIter, 2);
        clear lambdaShiftQRIter;

        % Compute approximate eigenvalues using Shifted QR Iteration w/ rel tol
        [lambdaShiftQRIterRel, timeShiftQRIterRel] = timeAlgorithm( ...
                                            algorithmShiftQRIterRel, T, n);
        timeShiftQRIterRel = timeShiftQRIterRel + timeTridiagonalize;
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
        n, minSpectralGap, errorNaiveBi, errorParBi, errorTreeBi, ...
        errorQRIter, errorShiftQRIter, errorShiftQRIterRel, errorJacobi, errorDivNConq, ...
        errorNaiveBi/normExact, errorParBi/normExact, errorTreeBi/normExact, ...
        errorQRIter/normExact, errorShiftQRIter/normExact, ...
        errorShiftQRIterRel/normExact, errorJacobi/normExact, errorDivNConq/normExact, ...
        timeNaiveBi, timeParBi, timeTreeBi, ...
        timeQRIter, timeShiftQRIter, timeShiftQRIterRel, timeJacobi, timeDivNConq, ...
        'VariableNames', ...
        {'n', 'MinSpectralGap', 'AbsError2Norm_NaiveBisect', 'AbsError2Norm_ParBisect', ...
         'AbsError2Norm_TreeBisect', 'AbsError2Norm_QRIter', ...
         'AbsError2Norm_ShiftQRIter', 'AbsError2Norm_ShiftQRIterRel', 'AbsError2Norm_Jacobi', ...
         'AbsError2Norm_DivNConq', ...
         'RelError2Norm_NaiveBisect', 'RelError2Norm_ParBisect', ...
         'RelError2Norm_TreeBisect', 'RelError2Norm_QRIter', ...
         'RelError2Norm_ShiftQRIter', 'RelError2Norm_ShiftQRIterRel', 'RelError2Norm_Jacobi', ...
         'RelError2Norm_DivNConq', ...
         'Time_NaiveBisect(in s)', 'Time_ParBisect(in s)', ...
         'Time_TreeBisect(in s)', 'Time_QRIter(in s)', ...
         'Time_ShiftQRIter(in s)', 'Time_ShiftQRIterRel(in s)', 'Time_Jacobi(in s)', 'Time_DivNConq(in s)'});

    % Clear intermediate vars for more space
    clear A Z T lambdaExact timeTridiagonalize normExact;

end


function [eigens, measuredTime] = timeAlgorithm(algorithm, Z, n)
    if n <= 256 % Use timeit for n <=256 since runtimes are relatively small
        measuredTime = timeit(@() algorithm(Z));
        eigens = algorithm(Z);
    else % Use tic/toc for n >= 256 since runtimes become large
        tic;
        eigens = algorithm(Z);
        measuredTime = toc;
    end
end