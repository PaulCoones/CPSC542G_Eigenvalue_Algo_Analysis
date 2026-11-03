%{
    compareAllStabilityMany.m Is a program which compares all relevant
    (eigenvalue) algorithms against one another on a non-perturbed and a
    perturbed tridiagonal matrix, by comparing their stability (change
    in output based on change in input), execution time, and errors
    within a csv table.

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
fileName = ['MassiveStabilityTests_' datestr(now, 'yyyymmdd_HHMMSS') '.csv'];

resultsFile = fullfile(outputDir, fileName);

% Test parameters (ensure length(nValues) = length(LValues))
n = 1024;
LValues = (n + 1) * 10.^(-3:3);
nValues = n * ones(size(LValues));

% Number of different perturbations to test for each (n, L) pair
numPerturbations = 5;

% Caclulate how many (n, L) pair tests we have to compute
numTests = length(nValues);

% Initalize the algorithms to compare using functions
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
    resultsRow = runBenchmark(n, L, numPerturbations, ...
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

% Delete the parallel pool finishing all parallel executions
delete(gcp('nocreate'));

fprintf('\n All %d tests completed. \n', numTests);
fprintf('All results saved to: \n %s \n', resultsFile);



function CSVresultsRow = runBenchmark(n, L, numPerturbations, ...
    algorithmDivNConq, algorithmTreeBi, algorithmParBi, ...
    algorithmNaiveBi, algorithmQRIter, algorithmShiftQRIter, ...
    algorithmShiftQRIterRel, algorithmJacobi)

    % Construct the (unperturbed) tridiagonal matrix
    T = triDiag(n, L);

    % Compute exact eigenvalues of the non-perturbed matrix
    h = L/(n+1);
    k = (1:n)';
    lambdaExactT = sort((4/h^2) * ...
        sin(k*pi/(2*(n+1))).^2);

    % Compute the minimum spectral gap of the non-perturbed matrix
    minSpectralGapT = min(diff(lambdaExactT));

    % Group algorithms for repeated perturbation tests
    algorithms = { ...
        algorithmNaiveBi, algorithmParBi, algorithmTreeBi, ...
        algorithmQRIter, algorithmShiftQRIter, ...
        algorithmShiftQRIterRel, algorithmJacobi, algorithmDivNConq};

    algorithmNames = { ...
        'Naive Bisection', 'Parallel Bisection', 'Tree Bisection', ...
        'QR Iteration', 'Shifted QR (Absolute)', ...
        'Shifted QR (Relative)', 'Jacobi', 'Divide and Conquer'};

    numAlgorithms = length(algorithms);

    % Initialize all errors for the current test
    errorDifference = zeros(1, numAlgorithms);
    errorPerturbed = zeros(1, numAlgorithms);
    errorNonPerturbed = zeros(1, numAlgorithms);

    % Initialize all timings for the current test
    timePerturbed = zeros(1, numAlgorithms);
    timeNonPerturbed = zeros(1, numAlgorithms);

    % Initialize the spectral gap and perturbation norm measurements
    minSpectralGapZ = 0;
    norm2DifferenceTZ = 0;

    % Initialize the sensitivity and perturbation response measurements
    SExact = 0;
    SAlgorithm = zeros(1, numAlgorithms);
    spectralChangeRatio = zeros(1, numAlgorithms);
    etaResponse = zeros(1, numAlgorithms);
    validRatioCount = 0;

    % Store the non-perturbed eigenvalues for each algorithm
    lambdaNonPerturbed = cell(1, numAlgorithms);

    % Compute approximate eigenvalues for the non-perturbed matrix
    % Only perform these computations once since T does not change
    for algorithmIndex = 1:numAlgorithms

        algorithm = algorithms{algorithmIndex};

        [lambdaT, timeNonPerturbed(algorithmIndex)] = ...
            timeAlgorithm(algorithm, T, n);

        % Incase eigenvalues not sorted, sort eigenvalues in ascending order
        lambdaT = sort(lambdaT(:));

        % Compute the error between estimated and exact non-perturbed eigenvalues
        errorNonPerturbed(algorithmIndex) = ...
            norm(lambdaExactT - lambdaT, 2);

        lambdaNonPerturbed{algorithmIndex} = lambdaT;

    end

    % Compute the perturbed tridiagonal matrix
    % Repeat the benchmark over multiple different perturbations
    for perturbationIndex = 1:numPerturbations

        fprintf('  Running perturbation %d of %d\n', ...
            perturbationIndex, numPerturbations);

        rng(7755 + perturbationIndex - 1);
        Z = T + perturbTriDiag(T);

        % Compute reference eigenvalues using MATLAB's eig solver
        lambdaExactZ = sort(eig(Z));

        % Compute the minimum spectral gap of the modified matrix
        minSpectralGapZ = minSpectralGapZ + ...
            min(diff(lambdaExactZ));

        % Compute the matrix 2-norm of the perturbation
        norm2DifferenceTZ = norm2DifferenceTZ + norm(T - Z, 2);

        % Compute the Frobenius norm of the perturbation
        perturbationNormFro = norm(Z - T, 'fro');

        if perturbationNormFro == 0
            error('The perturbation matrix is zero.');
        end

        % Compute the exact change in the eigenvalues
        exactSpectralDifference = lambdaExactZ - lambdaExactT;
        exactSpectralDifferenceNorm = norm(exactSpectralDifference, 2);

        % Compute the exact eigenvalue sensitivity
        currentSExact = exactSpectralDifferenceNorm ...
            / perturbationNormFro;

        SExact = SExact + currentSExact;

        % Check whether the exact spectral change is nonzero
        if exactSpectralDifferenceNorm > 0
            validRatioCount = validRatioCount + 1;
        else
            warning('Exact spectral change is zero. Ratio is undefined.');
        end

        % Compute approximate eigenvalues using each algorithm
        for algorithmIndex = 1:numAlgorithms

            algorithm = algorithms{algorithmIndex};
            lambdaT = lambdaNonPerturbed{algorithmIndex};

            % Compute approximate eigenvalues for the perturbed matrix
            [lambdaZ, measuredTime] = timeAlgorithm(algorithm, Z, n);

            % Incase eigenvalues not sorted, sort eigenvalues in ascending order
            lambdaZ = sort(lambdaZ(:));

            % Compute the difference between estimated perturbed and
            % estimated non-perturbed eigenvalues
            algorithmSpectralDifference = lambdaZ - lambdaT;

            currentErrorDifference = norm(algorithmSpectralDifference, 2);

            errorDifference(algorithmIndex) = ...
                errorDifference(algorithmIndex) + ...
                currentErrorDifference;

            % Compute the error between estimated and exact perturbed eigenvalues
            errorPerturbed(algorithmIndex) = ...
                errorPerturbed(algorithmIndex) + ...
                norm(lambdaExactZ - lambdaZ, 2);

            % Compute the algorithm's sensitivity to the perturbation
            currentSAlgorithm = currentErrorDifference ...
                / perturbationNormFro;

            SAlgorithm(algorithmIndex) = ...
                SAlgorithm(algorithmIndex) + currentSAlgorithm;

            % Compute the ratio between the algorithm's spectral change
            % and the exact spectral change
            if exactSpectralDifferenceNorm > 0
                currentSpectralChangeRatio = currentErrorDifference ...
                    / exactSpectralDifferenceNorm;

                spectralChangeRatio(algorithmIndex) = ...
                    spectralChangeRatio(algorithmIndex) + ...
                    currentSpectralChangeRatio;
            end

            % Compute the error in the algorithm's perturbation response
            currentEResponse = norm( ...
                algorithmSpectralDifference - exactSpectralDifference, 2);

            % Compute the normalized perturbation response error
            currentEtaResponse = currentEResponse / perturbationNormFro;

            etaResponse(algorithmIndex) = ...
                etaResponse(algorithmIndex) + currentEtaResponse;

            % Accumulate the timing for the perturbed matrix
            timePerturbed(algorithmIndex) = ...
                timePerturbed(algorithmIndex) + measuredTime;

            clear lambdaZ;

        end

        % Clear intermediate vars for more space
        clear Z lambdaExactZ;

    end

    % Compute averages over all perturbations
    minSpectralGapZ = minSpectralGapZ / numPerturbations;
    norm2DifferenceTZ = norm2DifferenceTZ / numPerturbations;

    errorDifference = errorDifference / numPerturbations;
    errorPerturbed = errorPerturbed / numPerturbations;
    timePerturbed = timePerturbed / numPerturbations;

    % Compute average sensitivity and perturbation response measurements
    SExact = SExact / numPerturbations;
    SAlgorithm = SAlgorithm / numPerturbations;
    etaResponse = etaResponse / numPerturbations;

    % Compute the average ratio over perturbations with nonzero spectral change
    if validRatioCount > 0
        spectralChangeRatio = spectralChangeRatio / validRatioCount;
    else
        spectralChangeRatio(:) = NaN;
    end

    % Tabulate the results for the currents test
    % Each row represents a measurement and each column an algorithm
    metricNames = { ...
        'n'; ...
        'L'; ...
        'NumPerturbations'; ...
        'MinSpectralGap_T'; ...
        'MeanMinSpectralGap_Z'; ...
        'MeanNorm2Difference_T_Z'; ...
        'MeanDifference2Norm_Eigenvalues'; ...
        'MeanS_exact'; ...
        'MeanS_algorithm'; ...
        'MeanSpectralChangeRatio'; ...
        'MeanEta_response'; ...
        'MeanError2Norm_Z'; ...
        'Error2Norm_T'; ...
        'MeanTime_Z(in s)'; ...
        'Time_T(in s)'};

    CSVresultsRow = table(string(metricNames), ...
        'VariableNames', {'Metric'});

    % Create one column per algorithm
    for algorithmIndex = 1:numAlgorithms

        CSVresultsRow.(algorithmNames{algorithmIndex}) = [ ...
            n; ...
            L; ...
            numPerturbations; ...
            minSpectralGapT; ...
            minSpectralGapZ; ...
            norm2DifferenceTZ; ...
            errorDifference(algorithmIndex); ...
            SExact; ...
            SAlgorithm(algorithmIndex); ...
            spectralChangeRatio(algorithmIndex); ...
            etaResponse(algorithmIndex); ...
            errorPerturbed(algorithmIndex); ...
            errorNonPerturbed(algorithmIndex); ...
            timePerturbed(algorithmIndex); ...
            timeNonPerturbed(algorithmIndex)];

    end

    % Clear intermediate vars for more space
    clear T lambdaExactT lambdaNonPerturbed k h;

end


function [eigens, measuredTime] = timeAlgorithm(algorithm, Z, n)
    if n > 4096 % Use timeit for n < 4096 since runtimes are relatively small
        measuredTime = timeit(@() algorithm(Z));
        eigens = algorithm(Z);
    else % Use tic/toc for n >= 4096 since runtimes become large
        tic;
        eigens = algorithm(Z);
        measuredTime = toc;
    end
end