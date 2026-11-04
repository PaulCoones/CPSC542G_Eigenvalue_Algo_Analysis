function deltaA = perturbTriDiag(A, seed, epsilon)
% x = triDiag(A) Returns a random symmetric tridiagonal matrix deltaA
% which has norm(deltaA) << norm(A), acting as a symmetric tridiagonal
% using the random seed provided.
% perturbation of A.
% Requires: A is a real square nxn matrix with n >= 1
%           epsilon is a real floating point number > 0 (ideally small)
%           seed is an integer

arguments
    A 
    seed = 0
    epsilon = 1e-3
end

% Initalizations
n = size(A, 1);
normA = norm(A, 2);
if seed ~= 0
    rng(seed);
end
deltaA = rand(n);

% Extract tridiagonal entries only
deltaA = triu(tril(deltaA, 1), -1);

% Enforce symmetry
deltaA = (deltaA + deltaA') / 2;

if normA == 0 % In-=case normA happens to be too small (below machine precs.)
    deltaA = zeros(n);
else
    deltaA = deltaA * (epsilon * normA / norm(deltaA, 2));
end

end