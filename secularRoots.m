function roots = secularRoots(D, beta, z, tol, k)

% MUST DEAL WITH DEFLATION BEFORE HEADING INTO A ROOT FINDING ALGORITHM!!!

% Extract and sort diagonal entries (FIX ME: DO I NEED TO SORT????)
d = sort(diag(D))

% Identify distinct poles
[poles, ~, groupInd] = unique(d) %unique(d, 'stable')?

% Count multiplicities by adding 1s according to the groupInds
multiplicity = accumarray(groupInd, 1) 

% Combine weights associated with repeated poles
weightSquared = accumarray(groupInd, z.^2)

% Retain groups with nonzero weights (according to 2-norm)
active = weightSquared > 0

% Known eigenvalues obtained from deflation
deflatedRoots = repelem( ...
    poles, multiplicity - double(active))

% Reduced d and z's
dReduced = poles(active);
zReduced = sqrt(weightSquared(active))

% d = diag(D); 
% 
% if beta == 0
%     warning("Vanishing rank-1 correction, " + ...
%         "eigenvalues already in the diagonal entries of D.")
%     roots = sort(d);
% end
% 
% 
% % Secular eqn and it's derivative
% f = @(lambda) 1 + beta * sum(z.^2 ./ (d - lambda));
% fp = @(lambda) beta * sum(z.^2 ./ (d - lambda).^2);
% 
% 
% w = sqrt(abs(beta)) * z; 




end