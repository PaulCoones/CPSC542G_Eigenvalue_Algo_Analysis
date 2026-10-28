function H = heatSource(x, y)
% s = sourceFcn(x, y) Returns the value of the source term, in this case
% representing two heating elements centered about (0.35, 0.6) and 
% (0.8, 0.25) both with radius 0.1.
% Requires: x and y are floating point numbers
% Runtime: O(1)

if ( sqrt((x - 0.35)^2 + (y - 0.6)^2) <= 0.1 || ...
        sqrt((x - 0.8)^2 + (y - 0.25)^2) <= 0.1 )
    H = 1;
else
    H = 0;
end

end
