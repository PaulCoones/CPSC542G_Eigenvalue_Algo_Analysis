%{ 
    Lap2D.m provides a standalone function for discretizing Poisson's equation
    subject to a source in 2D using finite differences, and thus producing a 
    matrix A and source vector b to be subsequently solved elsewhere.
%}

function [A, b] = Lap2D(m, L, sourceFcn)
% [A, b] = Lap2D(m) Returns an (m^2)x(m^2) matrix A and source vector b (mx1),
% corresponding to the discretization of poisson's equation with a source 
% term represented by the function sourceFcn using m grid-points in each dimension.
% Requires: m is a positive real integer
%           L is a positive real number
%           sourceFcn is a function handle
% Runtime: O(m^2)

dimA = m^2;  % Avoids recomputing m^2 multiple times
A = zeros(dimA, dimA);
b = zeros(dimA, 1);
h = L/(m+1); % Uniform grid spacing

for j = 1:m
    for i = 1:m
        k = i + (j - 1) * m; % Reaches m^2 for i = m and j = m as desired
        A(k, k) = 4;
        if i > 1
            A(k, k - 1) = -1; % left neighbor
        end
        if i < m
            A(k, k + 1) = -1; % right neighbor
        end
        if j > 1
            A(k, k - m) = -1; % lower neighbor (is m-values to left)
            % Note that we do not move down a row, since our 2D
            % grid is collapsed into 1-dimension. Rather we move across 1
            % grid length m in x  to represent a change in y-coordinate
        end
        if j < m
            A(k, k + m) = -1; % upper neighbor (is m-values to right)
        end
        x = i*h;
        y = j*h;
        b(k) = sourceFcn(x, y); % example source term
    end
end

A = A ./ h^2; % Can use /h^2 or ./h^2

end
