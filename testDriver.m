% FIX ME PAUL
% FIX ME PAUL

clear all;
close all;

% =====================================================
% INITIALIZATION
% =====================================================

nValues = [10; 25; 50; 60; 70; 80; 90; 100; ...
    250; 250; 325];

LValues = [ones(9,1); 250; 0.25];

numTests = length(nValues);
numMethods = 4;

% Eigenvalue algorithms
methods = {@eigenBisection, ...
           @eigenParBisection, ...
           @eigenTreeBisection, ...
           @eigenJacobi};

methodNames = {'NAIVE BISECTION', ...
               'PARALLEL BISECTION', ...
               'TREE BISECTION', ...
               'JACOBI'};

% Start parallel pool outside timing measurements
if isempty(gcp('nocreate'))
    parpool('Processes');
end

% Store matrices and eigenvalue vectors
Z = cell(numTests, 1);
lambda = cell(numTests, 1);
lambdaApprox = cell(numTests, numMethods);

% Initialize error norms for all methods
error2Norm = zeros(numTests, numMethods);
errorInfNorm = zeros(numTests, numMethods);
relative2Norm = zeros(numTests, numMethods);

% Initialize execution times
executionTime = zeros(numTests, numMethods);

% =====================================================
% COMPUTE EIGENVALUES AND ERRORS
% =====================================================

for i = 1:numTests

    n = nValues(i);
    L = LValues(i);

    % Construct tridiagonal matrix
    Z{i} = triDiag(n, L);

    % Initialize relevant eigenvalue vectors
    h = spacingVec(n, L);
    k = eigenStorage(n);

    % Compute exact eigenvalues
    lambda{i} = eigenTriExact(h, k, n);

    % Test each implementation independently
    for m = 1:numMethods

        % Time eigenvalue computation only
        tic;
        result = methods{m}(Z{i});
        executionTime(i,m) = toc;

        % Sort eigenvalues in ascending order
        lambdaApprox{i,m} = sort(result(:));

        % Compute error vector
        err = lambdaApprox{i,m} - lambda{i};

        % Compute error norms
        error2Norm(i,m) = norm(err, 2);
        errorInfNorm(i,m) = norm(err, inf);

        relative2Norm(i,m) = ...
            norm(err, 2) / norm(lambda{i}, 2);

    end
end

% =====================================================
% DISPLAY INDIVIDUAL RESULTS
% =====================================================

% Calculate percentage increases for each method
error2NormPctIncrease = NaN(numTests, numMethods);

for m = 1:numMethods
    if error2Norm(1,m) ~= 0
        error2NormPctIncrease(:,m) = ...
            ((error2Norm(:,m) - error2Norm(1,m)) ...
            / error2Norm(1,m)) * 100;
    end
end

% =====================================================
% NAIVE AND PARALLEL BISECTION RESULTS
% =====================================================

% Verify sequential and parallel implementations agree
for i = 1:numTests
    if ~isequaln(lambdaApprox{i,1}, lambdaApprox{i,2})
        warning('Sequential and parallel results differ for test %d.', i);
    end
end

% Both implementations use the same algorithm, so use
% naive bisection errors with separate execution times.
resultsBisection = table( ...
    nValues, LValues, ...
    error2Norm(:,1), ...
    errorInfNorm(:,1), ...
    relative2Norm(:,1), ...
    error2NormPctIncrease(:,1), ...
    executionTime(:,1), executionTime(:,2), ...
    'VariableNames', ...
    {'n', 'L', ...
     'Error_2Norm', 'Error_InfNorm', ...
     'Relative_2Norm', 'Error_2Norm_PctIncrease', ...
     'NaiveTime_s', 'ParallelTime_s'});

fprintf('\n=== NAIVE AND PARALLEL BISECTION RESULTS ===\n\n');

format short e
disp(resultsBisection);

% =====================================================
% TREE BISECTION RESULTS
% =====================================================

resultsTree = table( ...
    nValues, LValues, ...
    error2Norm(:,3), ...
    errorInfNorm(:,3), ...
    relative2Norm(:,3), ...
    error2NormPctIncrease(:,3), ...
    executionTime(:,3), ...
    'VariableNames', ...
    {'n', 'L', ...
     'Error_2Norm', 'Error_InfNorm', ...
     'Relative_2Norm', 'Error_2Norm_PctIncrease', ...
     'ExecutionTime_Seconds'});

fprintf('\n=== TREE BISECTION EIGENVALUE ERRORS ===\n\n');

disp(resultsTree);

% =====================================================
% JACOBI RESULTS
% =====================================================

resultsJacobi = table( ...
    nValues, LValues, ...
    error2Norm(:,4), ...
    errorInfNorm(:,4), ...
    relative2Norm(:,4), ...
    error2NormPctIncrease(:,4), ...
    executionTime(:,4), ...
    'VariableNames', ...
    {'n', 'L', ...
     'Error_2Norm', 'Error_InfNorm', ...
     'Relative_2Norm', 'Error_2Norm_PctIncrease', ...
     'ExecutionTime_Seconds'});

fprintf('\n=== JACOBI EIGENVALUE ERRORS ===\n\n');

disp(resultsJacobi);

% =====================================================
% EXECUTION TIME COMPARISON
% =====================================================

timeComparison = table( ...
    nValues, LValues, ...
    executionTime(:,1), executionTime(:,2), ...
    executionTime(:,3), executionTime(:,4), ...
    'VariableNames', ...
    {'n', 'L', ...
     'NaiveTime_s', 'ParallelTime_s', ...
     'TreeTime_s', 'JacobiTime_s'});

fprintf('\n=== EXECUTION TIME COMPARISON ===\n\n');

format short
disp(timeComparison);


% =====================================================
% EXECUTION TIME BENCHMARKS (25 RUNS)
% =====================================================

% Benchmark parameters
benchmarkN = [80; 250];
benchmarkL = [1; 250];

numRuns = 25;
numCases = length(benchmarkN);

% Store execution times:
% Dimensions: (run, algorithm, test case)
benchmarkTimes = zeros(numRuns, numMethods, numCases);

for i = 1:numCases

    n = benchmarkN(i);
    L = benchmarkL(i);

    % Construct test matrix
    A = triDiag(n, L);

    fprintf('\nBenchmarking n = %d, L = %g\n', n, L);

    % Warm up each algorithm before timing
    for m = 1:numMethods
        methods{m}(A);
    end

    % Perform 25 independent timing runs
    for r = 1:numRuns

        for m = 1:numMethods

            % Time algorithm execution only
            tStart = tic;
            methods{m}(A);
            benchmarkTimes(r,m,i) = toc(tStart);

        end

    end
end

% =====================================================
% BENCHMARK STATISTICS
% =====================================================

numRows = numCases * numMethods;

testN = zeros(numRows, 1);
testL = zeros(numRows, 1);
algorithm = strings(numRows, 1);

meanTime = zeros(numRows, 1);
medianTime = zeros(numRows, 1);
stdTime = zeros(numRows, 1);
minTime = zeros(numRows, 1);
maxTime = zeros(numRows, 1);

for i = 1:numCases
    for m = 1:numMethods

        row = (i-1)*numMethods + m;

        % Retrieve the 25 execution times
        times = benchmarkTimes(:,m,i);

        testN(row) = benchmarkN(i);
        testL(row) = benchmarkL(i);
        algorithm(row) = string(methodNames{m});

        % Compute statistics
        meanTime(row) = mean(times);
        medianTime(row) = median(times);
        stdTime(row) = std(times);
        minTime(row) = min(times);
        maxTime(row) = max(times);

    end
end

% =====================================================
% DISPLAY BENCHMARK RESULTS
% =====================================================

benchmarkResults = table( ...
    testN, testL, algorithm, ...
    meanTime, medianTime, stdTime, minTime, maxTime, ...
    'VariableNames', ...
    {'n', 'L', 'Algorithm', ...
     'Mean_s', 'Median_s', 'StdDev_s', ...
     'Min_s', 'Max_s'});

fprintf('\n=== EXECUTION TIME BENCHMARKS (25 RUNS) ===\n\n');

format short
disp(benchmarkResults);

% Store individual run times for further analysis
benchmark80 = benchmarkTimes(:,:,1);
benchmark250 = benchmarkTimes(:,:,2);


% =====================================================
% HELPER FUNCTIONS
% =====================================================

function h = spacingVec(n, L)
    h = L/(n+1);
end

function k = eigenStorage(n)
    k = (1:n)';
end

function lambda = eigenTriExact(h, k, n)
    lambda = (4/h^2) * sin(k*pi/(2*(n+1))).^2;
end