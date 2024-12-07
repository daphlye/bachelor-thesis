function [res] = NaNalyzer(input, max)
%NANALYZER Summary of this function goes here
% this function returns a 1 if less than 10% of the data in "input"-array is missing e.g.
% =NaN, else 0
    tot = numel(input);
    nans = sum(isnan(input));
    nanspertot = (nans / tot);
    res = nanspertot < max;
end
