function deltaA = perturbTriDiag(A)
% x = triDiag(A) Returns a random symmetric tridiagonal matrix deltaA
% which has norm(deltaA) << norm(A), acting as a symmetric tridiagonal
% perturbation of A.
% Requires: A is a real square nxn matrix with n >= 1

% Initalizations
n = size(A, 1);
normA = norm(A, 2);
rng(725911);
deltaA = rand(n);

% Extract tridiagonal entries only
deltaA = triu(tril(deltaA, 1), -1);

% Enforce symmetry
deltaA = (deltaA + deltaA') / 2;

% Ensure our randomized perturbation has a small relative norm
epsilon = 1e-3;

if normA == 0 % In-=case normA happens to be too small (below machine precs.)
    deltaA = zeros(n);
else
    deltaA = deltaA * (epsilon * normA / norm(deltaA, 2));
end

end