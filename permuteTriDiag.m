function deltaA = permuteTriDiag(A)
% x = triDiag(A) Returns a random symmetric tridiagonal matrix deltaA
% which has norm(deltaA) << norm(A).
% Requires: A is a real square nxn matrix with n >= 1

% Initalizations
n = size(A, 1);
normA = norm(A, 2);
deltaA = rand(n);

% Ensure the randomized permutation matrix has tiny norm
if normA < 0.01
    deltaA = deltaA * (0.01*normA);
else
    deltaA = deltaA / (100*normA);
end

deltaA = triu(tril(deltaA, 1), -1); % Extract diagonal entries only
deltaA = (deltaA + deltaA') / 2;    % Enforce symmetry

end