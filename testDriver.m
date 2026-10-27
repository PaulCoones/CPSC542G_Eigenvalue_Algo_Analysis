%{ 
    testDriver.m Is a program which compares all (eigenvalue)
    algorithms against one another, by tabulating their execution times
    and errors in their computed eigenvalues compared to the exact
    analytical eigenvalues.

    Hit run to run the program and display the table.
%}

% Clean workspace
clearvars;
close all;

% (n, L) Pairs to be tested
nValues = [15, 250, 250, 325, 500];
LValues = [1, 1, 251, 177, 501];
numTests = length(nValues);

% Change me if you want to change the number of workers
numWorkers = 8; 

% Eigenvalue algorithms
methodNames = {'Naive Bisection', ...
    'Parallel Naive Bisection', ...
    'Interval Tree Bisection', ...
    'Jacobi', ...
    'Divide and Conquer', ...
    'QR Algorithm'};

methods = {@eigenBisection, @eigenParBisection, ...
           @eigenTreeBisection, @eigenJacobi, ...
           @eigenDivNConq, @qreig};

numMethods = length(methods);

% Start the parallel pool
if isempty(gcp('nocreate'))
    parpool('Threads', numWorkers); % Can change to Processor(s) if desired
end

% Initialize variables for the table cells
numResults = numMethods * numTests;
algoNames = strings(numResults, 1);
execTime = zeros(numResults, 1);
norm2Err = zeros(numResults, 1);
nCol = zeros(numResults, 1); 
LCol = zeros(numResults, 1);

% Do work on each test pair (n, L)
for i = 1:numTests

    % Iterate over all testing pairs (n, L)
    n = nValues(i);
    L = LValues(i);

    % Construct tridiagonal matrix and compute all its n-known eigenvalues
    % Future considerations: Test other matrices other than triDiag
    T = triDiag(n, L);
    lambdaExact = exactTriDiagEigens(n/(L + 1), (1:n)', n); % k spans 1 -> n

    % For each fixed test pair (n, L), try each algorithm
    for j = 1:numMethods

        % Calculate the row index in our final table
        row = (i - 1) * numMethods + j;

        % Store current (n, L) pair and the current algo name
        nCol(row) = n;
        LCol(row) = L;
        algoNames(row) = string(methodNames{j});

        % Compute the eigenvalues by applying the current method to the
        % tridiagonal matrix T (sorted incase it is not already sorted)
        % and get the (2-norm) error against the exact solutions.
        lambdaApprox = sort(methods{j}(T));
        norm2Err(row) = norm(lambdaApprox - lambdaExact, 2);

        % Although the error is only from 1 computation, we do the timings
        % off various runs with timeit.
        execTime(row) = timeit(@() methods{j}(T));

    end
end

% Setup the table!
results = table(nCol, LCol, algoNames, execTime, norm2Err, ...
     'VariableNames', {'n', 'L', 'Algorithm', ...
     'Mean Time (s)', 'Mean 2-Norm Error'});

fprintf('\n=== EIGENVALUE ALGORITHM COMPARISON ===\n\n');

% Format (short scientific notation)
format short e

% Format n and L without scientific notation
results.n = compose("%.2f", nCol);
results.L = compose("%.2f", LCol);

% Display the table
disp(results);


% Returns the exact eigenvalue for the triDiagonal system provided the grid
% spacing h, the matrix size n, and the wave number k
function lambda = exactTriDiagEigens(h, k, n)
    lambda = (4/h^2) * sin(k*pi/(2*(n+1))).^2;
end