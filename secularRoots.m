function roots = secularRoots(D, beta, z, tol, k)

% Extract and sort diagonal entries (FIX ME: DO I NEED TO SORT????)
d = diag(D);
z = z(:)

if beta == 0
    roots = d;
    warning("Vanishing rank-1 correction, eigenvalues " + ...
            "already in the diagonal entries.")
    return;
end


% MUST DEAL WITH DEFLATION BEFORE HEADING INTO A ROOT FINDING ALGORITHM!!!
% Reference: SKIPPED in Trefethen & Bau's Numerical Linear Algebra Chapter 30
% So I referenced instead https://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapters5-6.pdf

% Identify UNIQUE poles of the "intended" secular function
[poles, ~, groupInd] = unique(d) %unique(d, 'stable')?

% Count multiplicities by adding 1s according to the groupInds
mult = accumarray(groupInd, 1) 

% Sum the values indices of z associated with repeated poles
weightSquared = accumarray(groupInd, z.^2)

% Retain groups with nonzero weights (according to 2-norm)
nonZero = weightSquared > 0

% Known eigenvalues obtained from deflation
% multiplicity - active term guarantees removal of non-repeated roots
deflatedRoots = repelem(poles, mult - nonZero)

% Adjusted d and z values with which we will define the secular equation
% and solve for its roots.
dAdj = poles(nonZero)
zAdj = sqrt(weightSquared(nonZero))

% Apply a mixed root finding algorithm combining Newton's and Bisection
% as Newton's cannot be guaranteed to converge, while Bisection can be.


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

% Try Newton's, if update diverges to infinity, resort to bisection
% Define the function we will solve
f = @(lambda) 1 + beta * sum(zAdj.^2 ./ (dAdj - lambda));
fp = @(lambda) beta * sum(zAdj.^2 ./ (dAdj - lambda).^2);

% Initalize
activeRoots = zeros(length(zAdj), 1);
numIntervals = length(zAdj);

% Call bisectionNewtons in parallel on all intervals
parfor intervalIndex = 1:numIntervals %  avoid recomputing length

    % Extract interval containing exactly one secular root
    left = searchIntervals(intervalIndex);
    right = searchIntervals(intervalIndex + 1);

    % Compute secular root using mixed Newton-Bisection algo
    activeRoots(intervalIndex) = bisectNewtons(f, fp, left, right, 1e-6, 100);

end

% Combine active secular roots with known deflated eigenvalues
roots = sort([activeRoots(:); deflatedRoots(:)]);

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
