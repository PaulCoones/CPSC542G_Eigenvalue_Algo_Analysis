%{
    compareTriDiagVaryL.m Is a program which compares two
    (eigenvalue) algorithms on the tridiagonal matrix constructed from
    finite differences. n is fixed at 128 and L (length of the rod) is 
    varied to control the spectral gap.

    Hit run to run the program and save the results to a CSV table.
    Change the Test parameters and/or algorithm parameters to compare
    different L values and/or algorithms.

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
fileName = ['compareTwoTriDiagVaryL_' datestr(now, 'yyyymmdd_HHMMSS') '.csv'];

resultsFile = fullfile(outputDir, fileName);

% Test parameters (ensure length(nValues) = length(LValues))
LValues = (128 + 1) * 10.^(-3:3);    % L = 0.129 to 129000 (h = 1e-3 to 1e3)
nValues = 128 * ones(size(LValues));  % Keep n fixed at 128 for every test

% Caclulate how many (n, L) pair tests we have to compute
numTests = length(nValues);

% Initalize the two algorithms to compare with function handles
algorithmX = @eigenJacobi;        % Algorithm X
algorithmY = @eigenDivNConqNew;   % Algorithm Y

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
    resultsRow = runBenchmark(n, L, algorithmX, algorithmY);

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


function CSVresultsRow = runBenchmark(n, L, algorithmX, algorithmY)

    % Construct the tridiagonal matrix
    Z = triDiag(n, L);

    % Compute all exact eigenvalues from the analytical solution
    h = L/(n+1);
    k = (1:n)';
    lambdaExact = (4/h^2) * ...
        sin(k*pi/(2*(n+1))).^2;

    % Spectral gap between the two smallest exact eigenvalues
    spectralGap = lambdaExact(2) - lambdaExact(1);

    % Compute approximate eigenvalues using Algorithm X
    [lambdaX, timeX] = timeAlgorithm(algorithmX, Z, n);

    % Incase eigenvalues not sorted, sort eigenvalues in ascending order
    lambdaX = sort(lambdaX(:));
    errorX = norm(lambdaExact - lambdaX, 2);        % Absolute error
    relErrorX = errorX / norm(lambdaExact, 2);      % Relative error
    clear lambdaX; % Clear memory


    % Compute approximate eigenvalues using Algorithm Y
    [lambdaY, timeY] = timeAlgorithm(algorithmY, Z, n);
    lambdaY = sort(lambdaY(:));
    errorY = norm(lambdaExact - lambdaY, 2);
    relErrorY = errorY / norm(lambdaExact, 2);
    clear lambdaY;

    % Tabulate the results for the currents test
    CSVresultsRow = table( ...
        n, L, spectralGap, errorX, errorY, relErrorX, relErrorY, ...
        timeX, timeY, ...
        'VariableNames', ...
        {'n', 'L', 'SpectralGap', ...
         'Error2Norm_AlgoX', 'Error2Norm_AlgoY', ...
         'RelError2Norm_AlgoX', 'RelError2Norm_AlgoY', ...
         'Time_AlgoX(in s)', 'Time_AlgoY(in s)'});

    % Clear intermediate vars for more space
    clear Z lambdaExact k h;

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