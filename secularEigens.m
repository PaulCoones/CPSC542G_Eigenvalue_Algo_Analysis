function [roots, eigenVectors] = secularEigens(D, beta, z, tol, max_iters)

% Extract diagonal entries for ease of readability
d = diag(D);
n = length(D);

% Ensure we need to actually do work
if beta == 0
    warning("Vanishing rank-1 correction, eigenvalues " + ...
        "already in the diagonal entries.")

    % Return eigenvalues and eigenvectors in increasing order
    [roots, newOrdering] = sort(d);
    eigenVectors = eye(n, n);
    eigenVectors = eigenVectors(:, newOrdering); % follows the sorting above
    return;
end
    
% MUST DEAL WITH DEFLATION BEFORE HEADING INTO A ROOT FINDING ALGORITHM!!!
% Reference: SKIPPED in Trefethen & Bau's Numerical Linear Algebra Chapter 30
% So I referenced instead https://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapters5-6.pdf

% Identify UNIQUE poles of the "intended" secular function
[poles, ~, groupInd] = unique(d);
numPoles = length(poles); % Record the number of poles for later use

% Count multiplicities by adding 1s according to the groupInds
mult = accumarray(groupInd, 1);

% Sum the values indices of z associated with repeated poles
weightSquared = accumarray(groupInd, z.^2);

% Retain groups with nonzero weights (according to 2-norm)
nonZero = weightSquared > 0;

% Known eigenvalues obtained from deflation
% multiplicity - active term guarantees removal of non-repeated roots
deflatedRoots = repelem(poles, mult - nonZero);

% Adjusted d and z values with which we will define the secular equation
% and solve for its roots.
dAdj = poles(nonZero);
zAdj = sqrt(weightSquared(nonZero));

% Construct a cumulative Givens rotations matrix to obtain the 
% eigenvectors associated with the deflated eigenvalues.
G = eye(n, n);

% Construct a cumulative (hence temporary) rotated version of z
zTemp = z;

% Keep track of the indics of the deflated eigenvectors
deflatedInds = zeros(n - length(dAdj), 1);

% Keep track of the deflated eigenvectors themseleves
deflatedEigens = zeros(n, n - length(dAdj));

% Keep track of the number of deflated eigenvectors dealt with so far
curDeflateCount = 0;

% Do work on all the poles
for pole = 1:numPoles

    % Get all indices associated with the current pole
    curPoleInds = find(pole == groupInd); 
    numPoleInds = length(curPoleInds); % How many repeats

    % Check is all components of z corresponding to the pole sum to zero
    % If they do, then every corresponding standard basis vector e_i 
    % is an eigenvector of D + beta*z*z' = D on that subspace, 
    % and we can stop looking at this pole immediately and move on
    if ~nonZero(pole) 
        deflatedInds(curDeflateCount + (1:numPoleInds)) = curPoleInds;
        curDeflateCount = curDeflateCount + numPoleInds;

        continue; % Move onto next iteration as we are done here
    end

    % Set a pivot around which we will rotate
    curPivotInd = curPoleInds(1);
    
    % Do extra work on repeated poles with Givens Rotations
    for k = 2:numPoleInds

        % Adjust the current pivot according to the new value given by temp
        curPivotVal = zTemp(curPivotInd);
            
        % Choose the component to zero out
        indRotatedOut = curPoleInds(k);
        valRotatedOut = zTemp(indRotatedOut);

        % Find the Givens Rotation angles
        newPivotVal = sqrt(curPivotVal*curPivotVal + ...
                           valRotatedOut*valRotatedOut); 
        % WHAT IF ZERO? IS IT POSSIBLEz?????

        c = curPivotVal / newPivotVal;
        s = valRotatedOut / newPivotVal; 

        % Apply the current Givens rotation to zero out the corresponding
        % component in the cumulative Givens rotation matrix G.
        G(:, [curPivotInd, indRotatedOut]) = G(:, [curPivotInd, indRotatedOut]) * ...
                                             [c, -s; s, c];

        % Adjust z so that it now has a zeroed out component toRotateOut
        zTemp([curPivotInd, indRotatedOut]) = [newPivotVal; 0];
    end

    % Calculate how many deflations we've done (all but the pivot)
    numDeflations = length(curPoleInds) - 1;

    if numDeflations > 0
        % Keep track of which indices we have deflated to later determine
        % eigenvectors
        deflatedInds(curDeflateCount + (1:numDeflations)) = curPoleInds(2:end);

        % Adjust number of deflations completed in total so far
        curDeflateCount = curDeflateCount + numDeflations; 
    end

end

% Get the eigenvectors, equivalent to G * I(:, deflatedInds), 
% with I = eye(n,n).
deflatedEigens = G(:, deflatedInds);

% Lowest possible value of any eigenvalue by Gershgorin's
% Recalling M[i,i] = D[i,i] + beta*z[i]*z[i]
% While     M[i,j] = beta*z[i]*z[j]
% Note: sum(abs(zAdj)) - abs(zAdj) produces sum_{i \neq j} zAdj[i]
gershCentres = dAdj + beta*zAdj.^2;
gershRadii = abs(beta)*abs(zAdj) .* (sum(abs(zAdj)) - abs(zAdj));

% Create m search intervals each interval containing ONLY ONE root
if beta > 0
    % Interlacing follows:
    % dAdj[1] < lambda_1 < dAdj[2] < ... < dAdj[p] < lambda_p < upperBound
    upperBound = max(gershCentres + gershRadii);
    searchIntervals = [dAdj; upperBound];
else
    % Interlacing follows:
    % lowerBound < lambda_1 < dAdj[1] < lambda_2 < ... < lambda_p < dAdj[p]
    lowerBound = min(gershCentres - gershRadii);
    searchIntervals = [lowerBound; dAdj];
end

% Define the "deflated" secular equation which we will solve 
% and its deriv as required by Newton's method.
f = @(lambda) 1 + beta * sum(zAdj.^2 ./ (dAdj - lambda));
fp = @(lambda) beta * sum(zAdj.^2 ./ (dAdj - lambda).^2);

% Initalize
nonDeflatedRoots = zeros(length(zAdj), 1);
numIntervals = length(zAdj);

% Call a root finding algorithm in parallel on all intervals
parfor intervalIndex = 1:numIntervals

    % Extract intervals, each interval which contains ONLY ONE root
    left = searchIntervals(intervalIndex);
    right = searchIntervals(intervalIndex + 1);

    % Compute root using a root finding algoritm
    nonDeflatedRoots(intervalIndex) = bisectNewtons(f, fp, left, right, ...
                                                    tol, max_iters);

end

% Explicitly construct the eigenvectors for the nonDeflatedRoots
% In particular, the eigenvectors of the matrix 
% MAdj = diag(dAdj) + beta * (zAdj*zAdj');
% Why can I use numIntervals? Well every interval contains exactly
% one eigenvalue, so the number of eigenvectors/eigenvalues = number of
% intervals.
nonDeflatedEigenVecs = zeros(numIntervals, numIntervals);

for rootInd = 1:numIntervals
    eigenvector = zAdj ./ (dAdj - nonDeflatedRoots(rootInd)); % elem-wise
    nonDeflatedEigenVecs(:, rootInd) = eigenvector / norm(eigenvector);
end

% Obtain the non-deflated coordinate indices
nonDeflatedInds = setdiff((1:n)', deflatedInds);

% MUST REORDER TO MATCH INCREASING POLE ORDER OF dAdj var
[~, newOrdering] = sort(groupInd(nonDeflatedInds)); % Only need the indices, not the sorted version
nonDeflatedInds = nonDeflatedInds(newOrdering);

% Map to the the original coordinates
V = G(:, nonDeflatedInds) * diag(sign(zTemp(nonDeflatedInds)));

% Transform the reduced eigenvectors
nonDeflatedEigenVecs = V * nonDeflatedEigenVecs;

% Combine the deflated and non deflated roots
[roots, newOrdering] = sort([nonDeflatedRoots(:); deflatedRoots(:)]);

% Combine deflated and non deflated eigenvectors
eigenVectors = [nonDeflatedEigenVecs, deflatedEigens];

% Reorder eigenvectors to match the ordering of the eigenvalues
% (increasingly sorted by sort(...)).
eigenVectors = eigenVectors(:, newOrdering);

end



function root = bisectNewtons(f, fp, L, U, tol, maxIters)

% Choose an initial guess in the middle of the interval
guess = L + (U - L)/2; 

% Check for intervals that are too small
if guess == L || guess == U
    if (U - L) <= tol
        root = guess;
        return;
    else 
        error(['Machine precision limitations is causing evaluation at' ...
            'the poles, which we cannot do.']);
    end
end

% After guaranteeing we are not at a pole, evaluate the secular eqn at our
% initial guess
f_iter = f(guess);

% Quick exit if our initial guess is sufficient
if f_iter == 0
    root = guess;
    return;
end

direction = sign(fp(guess));

if ~isfinite(direction) || direction == 0
    error('Cannot determine secular function monotonicity.');
end

k = 0;
while (k < maxIters)

    % Increment
    k = k + 1; 

    % If not yet finished.
    % Update the interval no matter if the next guess was obtained
    % from Newton's method or Bisection.
    % We cannot check opposite signed endpoints like in Bisection,
    % as we cannot evaluate the secular eqn at the poles.
    % Instead, we use the fact that the secular equation is STRICTLY
    % MONOTONIC between poles, then for example, if
    % f_iter < 0 and fp_iter > 0, then we can replace our negative
    % left endpoint L with f_iter, which is also negative.
    % Similar results follow for the 3 other cases.
    if (direction * sign(f_iter)) < 0
        L = guess;
    else
        U = guess;
    end

    % Check root tolerance by seeing how small our interval is
    if (U - L) <= tol
        root = L + (U - L)/2;
        return;
    end

    % Update f'
    fp_iter = fp(guess);

    % Perform a step of Newton's if it is safe to do so
    if (isfinite(fp_iter) && fp_iter ~= 0 && isfinite(f_iter))
        nextGuess = guess - (f_iter / fp_iter);
    else
        nextGuess = NaN; 
    end

    % If the step takes us to infinity, or is less/greater than the
    % left/right endpoints respectively, resort to bisection instead
    if (~isfinite(nextGuess) || nextGuess <= L || ...
        nextGuess >= U || nextGuess == guess) 
        
        % Bisect the interval
        nextGuess = L + (U - L)/2;
        
    end

    % Check if machine precision has prevented the guess from moving
    % forward at all, causing potential evaluations at the poles of the 
    % secular equation
    if nextGuess == L || nextGuess == U
        error('Machine precision causing divergence.');
    end

    guess = nextGuess;

    if isnan(guess)
        error('Guess is NaN.') % Should never reach this
    end

    % Update f
    f_iter = f(guess);

    if f_iter == 0
        root = guess;
        return;
    end

end

warning('Maximum iterations hit before convergence criterion met.');

end


% % CONSOLE TEST FOR THE PROGRAM
% % 12x12 diagonal matrix with repeated eigenvalues
% D = diag([8; -3; 2; 8; 0; 15; -3; 5; 0; 22; 11; 5]);
% 
% % Rank-one perturbation parameters
% beta = 1.75;
% z = [1.2; -0.9; 0; 0.4; 1.5; -1.1; ...
%     2.0; 0; -0.3; 0.8; -1.4; 0.6];
% 
% % Compute eigenvalues using your secular equation solver
% computedRoots = secularRoots(D, beta, z, 1e-9, 300);
% 
% % Compute reference eigenvalues using MATLAB
% referenceRoots = sort(eig(D + beta*(z*z')));
% 
% % Compare the results
% disp(table(computedRoots, referenceRoots, ...
%     abs(computedRoots - referenceRoots), ...
%     'VariableNames', {'Computed', 'MATLAB', 'AbsoluteError'}));
% 
% % Maximum absolute error
% fprintf('Maximum absolute error: %.12e\n', ...
%     max(abs(computedRoots - referenceRoots)));